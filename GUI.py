import os
import sys
import csv
import re
import queue
import threading
from datetime import datetime
from pathlib import Path
from tkinter import Tk, ttk, StringVar, END, DISABLED, N, S, E, W
from tkinter import messagebox, scrolledtext

# ====== User-tweakables ======
BAUD = 115200
READ_TIMEOUT = 0.1  # seconds
GUI_POLL_MS = 50    # ms
SESSIONS_DIRNAME = "PD_Sessions"  # Desktop folder for CSVs

# New: DAC clamp limits (edit if needed)
DAC_MIN_V = 0.0
DAC_MAX_V = 4.9      # keep a little headroom below 5V
DAC_OFFSET_V = 0.01   # we will set DAC to (A0 - 0.5 V)

# ====== Known names (helpful labels) ======
I2C_NAMES = {
    "0x40": "TLC59108 LED Driver",
    "0x4F": "MAX30205 Temp",
    "0x4C": "DAC7571",
    "0x4D": "DAC7571",
    "0x4E": "DAC7571",
    "0x4F": "DAC7571",
}

# ====== Parsing patterns (robust) ======
PATTERNS = {
    "autorun_begin": re.compile(r"^===\s*Full\s*Cycle\s*Run", re.I),
    "autorun_end":   re.compile(r"^===\s*Autorun\s*complete\s*===\s*$", re.I),
    "channel":       re.compile(r"^===\s*Channel\s+(\d+)\s*===\s*$", re.I),

    "temp":          re.compile(r"^TEMP:\s*([\-NaN\d\.]+)", re.I),
    "temp_on":       re.compile(r"^TEMP\s*\(on\):\s*([\-NaN\d\.]+)", re.I),

    "a0":            re.compile(r"^A0:\s*([\-NaN\d\.]+)", re.I),
    "a1":            re.compile(r"^A1:\s*([\-NaN\d\.]+)", re.I),
    "a0_dark":       re.compile(r"^A0\s*\(dark\):\s*([\-NaN\d\.]+)", re.I),
    "a0_on":         re.compile(r"^A0\s*\(on\):\s*([\-NaN\d\.]+)", re.I),
    "a0_post":       re.compile(r"^A0\s*\(post\-off\):\s*([\-NaN\d\.]+)", re.I),
    
    "a1_dark":       re.compile(r"^A1\s*\(dark\):\s*([\-NaN\d\.]+)", re.I),
    "a1_on":         re.compile(r"^A1\s*\(on\):\s*([\-NaN\d\.]+)", re.I),
    "a1_post":       re.compile(r"^A1\s*\(post\-off\):\s*([\-NaN\d\.]+)", re.I),


    "dv":            re.compile(r"^ΔV\s*=\s*([\-NaN\d\.]+)", re.I),

    "i2c_found":     re.compile(r"^(.*)\s:\s0x([0-9A-Fa-f]{2})$"),
    "i2c_found_flex": re.compile(r"(?:found.*at\s*)?0x([0-9A-Fa-f]{2})", re.I),

    "dac_ack":       re.compile(r"^DAC(?:7571)?\s*set\s*to\s*([\-NaN\d\.]+)\s*V", re.I),
    "dac_raw_ack":   re.compile(r"^DAC\s*raw\s*0x([0-9A-Fa-f]{2})\s+(\d+)", re.I),
}

CSV_HEADERS = ["timestamp_iso", "command_sent", "event", "channel", "value", "units", "raw_line"]

# ====== Serial tools ======
def list_serial_ports():
    try:
        import serial.tools.list_ports
        return [p.device for p in serial.tools.list_ports.comports()]
    except Exception:
        return []

class SerialWorker(threading.Thread):
    def __init__(self, port, baud, line_queue, on_closed):
        super().__init__(daemon=True)
        self.port = port
        self.baud = baud
        self.queue = line_queue
        self.on_closed = on_closed
        self._stop = threading.Event()
        self.ser = None

    def run(self):
        try:
            import serial
            self.ser = serial.Serial(self.port, self.baud, timeout=READ_TIMEOUT)
            while not self._stop.is_set():
                try:
                    line = self.ser.readline()
                    if line:
                        text = line.decode(errors="replace").strip()
                        self.queue.put(text)
                except Exception:
                    break
        except Exception as e:
            self.queue.put(f"[ERROR] Could not open {self.port}: {e}")
        finally:
            if self.ser:
                try: self.ser.close()
                except Exception: pass
            self.on_closed()

    def send(self, s: str):
        if self.ser and self.ser.is_open:
            self.ser.write(s.encode())

    def stop(self):
        self._stop.set()

# ====== CSV Logger ======
class CsvLogger:
    def __init__(self):
        desktop = Path.home() / "Desktop"
        target_dir = desktop / SESSIONS_DIRNAME
        target_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.path = target_dir / f"pd_session_{ts}.csv"
        self._file = open(self.path, "w", newline="", encoding="utf-8")
        self.writer = csv.writer(self._file)
        self.writer.writerow(CSV_HEADERS)
        self._file.flush()
        self.last_command = ""

    def set_last_command(self, cmd: str):
        self.last_command = cmd

    def log(self, event, channel, value, units, raw_line):
        self.writer.writerow([
            datetime.now().isoformat(timespec="seconds"),
            self.last_command,
            event,
            channel if channel is not None else "",
            value if value is not None else "",
            units or "",
            raw_line
        ])
        self._file.flush()

    def close(self):
        try: self._file.close()
        except Exception: pass

# ====== GUI App ======
class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Photodiode Board Controller — Serial GUI")
        self.root.geometry("1040x680")

        # State
        self.serial_worker = None
        self.line_queue = queue.Queue()
        self.csv = CsvLogger()
        self.current_channel = None

        # Autorun aggregation state
        self.autorun_active = False
        self._summary_emitted = False   # one-shot guard for autorun summary

        self.autorun_dv = {}        # channel -> ΔV (V)
        
                # A0 (pre-diff) deltas are filled from "ΔV = ..." line the firmware prints
        self.autorun_dv_a0 = {}     # channel -> ΔV (V) for A0 (pre-diff)
        # A1 (post-diff) deltas we compute as A1_on - A1_dark
        self.autorun_dv_a1 = {}     # channel -> ΔV (V) for A1 (post-diff)
        self._a1_dark = {}          # channel -> last A1_dark

        # Reference margin (V) the user can edit (used for DAC = A0_dark - margin)
        self.ref_margin_var = StringVar(value="0.500")

        
        # New: track DAC actions per channel per autorun, and last A0 seen
        self._dac_done_for_channel = set()
        self._last_a0_seen = {}     # channel -> float

        # Devices seen (addr -> name)
        self.devices_seen = {}

        # ===== Top bar =====
        top = ttk.Frame(root, padding=8)
        top.grid(row=0, column=0, sticky=E+W)
        ttk.Label(top, text="Serial Port:").grid(row=0, column=0, padx=(0,6))
        self.port_var = StringVar()
        self.port_combo = ttk.Combobox(top, textvariable=self.port_var, width=30, state="readonly")
        self.port_combo.grid(row=0, column=1, padx=(0,6))
        self._refresh_ports()
        ttk.Button(top, text="↻ Refresh", command=self._refresh_ports).grid(row=0, column=2, padx=(0,6))
        self.btn_connect = ttk.Button(top, text="Connect", command=self.connect)
        self.btn_connect.grid(row=0, column=3, padx=(0,6))
        self.btn_disconnect = ttk.Button(top, text="Disconnect", command=self.disconnect, state=DISABLED)
        self.btn_disconnect.grid(row=0, column=4, padx=(0,6))
        self.btn_open_csv = ttk.Button(top, text="Open CSV Folder", command=self.open_csv_folder)
        self.btn_open_csv.grid(row=0, column=5)
        self._append_notice = f"[CSV] Logging to: {self.csv.path}\n"

        # ===== Commands row =====
        btns = ttk.LabelFrame(root, text="Commands", padding=8)
        btns.grid(row=1, column=0, sticky=E+W, padx=8, pady=(0,8))
        self._mk_button(btns, "Full Cycle Run (r)", lambda: self.send_cmd("r"), 0, 0)
        self._mk_button(btns, "Temp (t)",          lambda: self.send_cmd("t"), 0, 1)
        self._mk_button(btns, "A0 (a0)",           lambda: self.send_cmd("a0"), 0, 2)
        self._mk_button(btns, "A1 (a1)",           lambda: self.send_cmd("a1"), 0, 3)
        self._mk_button(btns, "I2C Scan",          self.scan_i2c,               0, 4)
        self._mk_button(btns, "Help/menu",         lambda: self.send_cmd("help"), 0, 5)

        led_on = ttk.LabelFrame(btns, text="LED ON @255")
        led_on.grid(row=1, column=0, columnspan=3, sticky=E+W, pady=(8,0))
        for ch in range(4):
            self._mk_button(led_on, f"LED {ch} ON",  lambda c=ch: self.send_cmd(f"{c} on"), 0, ch)

        led_off = ttk.LabelFrame(btns, text="LED OFF")
        led_off.grid(row=1, column=3, columnspan=3, sticky=E+W, pady=(8,0))
        for ch in range(4):
            self._mk_button(led_off, f"LED {ch} OFF", lambda c=ch: self.send_cmd(f"{c} off"), 0, ch)

        # ===== DAC panel =====
        dac = ttk.LabelFrame(root, text="DAC7571 — Set Output", padding=8)
        dac.grid(row=2, column=0, sticky=E+W, padx=8, pady=(0,8))

        ttk.Label(dac, text="Address:").grid(row=0, column=0, sticky=E, padx=4, pady=2)
        self.dac_addr_var = StringVar(value="")
        self.dac_addr_combo = ttk.Combobox(dac, textvariable=self.dac_addr_var, width=10, state="normal")
        self.dac_addr_combo["values"] = ["", "0x4C", "0x4D", "0x4E", "0x4F"]
        self.dac_addr_combo.grid(row=0, column=1, sticky=W, padx=4, pady=2)

        ttk.Label(dac, text="Voltage (V):").grid(row=0, column=2, sticky=E, padx=4, pady=2)
        self.dac_volt_var = StringVar(value="2.500")
        volt_entry = ttk.Entry(dac, textvariable=self.dac_volt_var, width=10)
        volt_entry.grid(row=0, column=3, sticky=W, padx=4, pady=2)
        
        ttk.Label(dac, text="Ref margin (V):").grid(row=0, column=5, sticky=E, padx=4, pady=2)
        self.ref_margin_var = getattr(self, "ref_margin_var", StringVar(value="0.500"))
        ref_entry = ttk.Entry(dac, textvariable=self.ref_margin_var, width=10)
        ref_entry.grid(row=0, column=6, sticky=W, padx=4, pady=2)


        self.btn_dac_set = ttk.Button(dac, text="Set DAC (Volts)", command=self._send_dac_volts)
        self.btn_dac_set.grid(row=0, column=4, sticky=W, padx=6, pady=2)

        ttk.Label(dac, text="Raw 12-bit code:").grid(row=1, column=0, sticky=E, padx=4, pady=2)
        self.dac_raw_code_var = StringVar(value="")
        raw_entry = ttk.Entry(dac, textvariable=self.dac_raw_code_var, width=10)
        raw_entry.grid(row=1, column=1, sticky=W, padx=4, pady=2)

        self.btn_dac_raw = ttk.Button(dac, text="Send Raw Code", command=self._send_dac_raw)
        self.btn_dac_raw.grid(row=1, column=2, sticky=W, padx=6, pady=2)

        ttk.Label(dac, text="Detected DAC(s):").grid(row=1, column=3, sticky=E, padx=4, pady=2)
        self.detected_dac_var = StringVar(value="(none)")
        ttk.Label(dac, textvariable=self.detected_dac_var).grid(row=1, column=4, sticky=W, padx=4, pady=2)

        # ===== Devices panel =====
        devs = ttk.LabelFrame(root, text="I²C Devices (from last scan)", padding=8)
        devs.grid(row=3, column=0, sticky=E+W, padx=8, pady=(0,8))
        cols = ("addr", "name")
        self.devs_tree = ttk.Treeview(devs, columns=cols, show="headings", height=6)
        self.devs_tree.heading("addr", text="Address")
        self.devs_tree.heading("name", text="Name")
        self.devs_tree.column("addr", width=90, anchor="center")
        self.devs_tree.column("name", width=300, anchor="w")
        self.devs_tree.grid(row=0, column=0, sticky=E+W)
        devs.grid_columnconfigure(0, weight=1)

        # ===== Console =====
        self.console = scrolledtext.ScrolledText(root, wrap="word", height=16)
        self.console.grid(row=4, column=0, sticky=N+S+E+W, padx=8, pady=(0,8))
        self._append_console(self._append_notice)

        # Layout stretch
        root.grid_rowconfigure(4, weight=1)
        root.grid_columnconfigure(0, weight=1)

        # Poll serial
        self.root.after(GUI_POLL_MS, self._poll_serial)

        # Close handler
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

    # --- UI helpers ---
    def _mk_button(self, parent, text, cmd, r, c):
        b = ttk.Button(parent, text=text, command=cmd)
        b.grid(row=r, column=c, padx=4, pady=4, sticky=E+W)
        return b

    def _refresh_ports(self):
        ports = list_serial_ports()
        self.port_combo["values"] = ports
        if ports and (self.port_var.get() == "" or self.port_var.get() not in ports):
            self.port_var.set(ports[0])

    def connect(self):
        port = self.port_var.get()
        if not port:
            messagebox.showwarning("Select Port", "Please choose a serial port first.")
            return
        if self.serial_worker:
            messagebox.showinfo("Already Connected", "Serial is already connected.")
            return

        def on_closed():
            self.root.after(0, self._on_serial_closed)

        self.line_queue.queue.clear()
        self.serial_worker = SerialWorker(port, BAUD, self.line_queue, on_closed)
        self.serial_worker.start()
        self._append_console(f"[INFO] Connecting to {port} @ {BAUD}...\n")
        self.btn_connect.config(state=DISABLED)
        self.btn_disconnect.config(state="normal")

    def disconnect(self):
        if self.serial_worker:
            self.serial_worker.stop()
            self.serial_worker = None

    def _on_serial_closed(self):
        self._append_console("[INFO] Serial closed.\n")
        self.btn_connect.config(state="normal")
        self.btn_disconnect.config(state=DISABLED)

    def send_cmd(self, cmd: str):
        if cmd.strip().lower() == "r":
            self.autorun_active = True
            self.autorun_dv.clear()
            self.autorun_dv_a0.clear()
            self.autorun_dv_a1.clear()
            self._a1_dark.clear()
            self._dac_done_for_channel.clear()
            self._last_a0_seen.clear()


        if not self.serial_worker:
            messagebox.showwarning("Not Connected", "Connect to a serial port first.")
            return
        line = cmd.strip() + "\n"
        try:
            self.serial_worker.send(line)
            self.csv.set_last_command(cmd.strip())
            self._append_console(f"> {cmd.strip()}\n")
        except Exception as e:
            self._append_console(f"[ERROR] send failed: {e}\n")

    def scan_i2c(self):
        self.devices_seen.clear()
        self._refresh_devices_panel()
        self.detected_dac_var.set("(scanning...)")
        self.send_cmd("scan")

    def _poll_serial(self):
        try:
            while True:
                line = self.line_queue.get_nowait()
                self._handle_line(line)
        except queue.Empty:
            pass
        self.root.after(GUI_POLL_MS, self._poll_serial)

    def _append_console(self, text: str):
        self.console.insert(END, text)
        self.console.see(END)

    def open_csv_folder(self):
        folder = Path(self.csv.path).parent
        try:
            os.startfile(folder)  # Windows
        except Exception:
            try:
                if sys.platform == "darwin":
                    os.system(f"open '{folder}'")
                else:
                    os.system(f"xdg-open '{folder}'")
            except Exception:
                messagebox.showinfo("CSV Folder", f"Folder:\n{folder}")

    def on_close(self):
        self.disconnect()
        self.csv.close()
        self.root.destroy()

    # --- DAC controls ---
    def _send_dac_volts(self):
        addr = self.dac_addr_var.get().strip()
        vtxt = self.dac_volt_var.get().strip()
        if not vtxt:
            messagebox.showwarning("DAC Voltage", "Enter a target voltage, e.g., 2.500")
            return
        try:
            float(vtxt)
        except ValueError:
            messagebox.showwarning("DAC Voltage", "Voltage must be a number.")
            return

        if addr:
            self._raw_send(f"dac {addr} {vtxt}")
        else:
            self._raw_send(f"dac {vtxt}")

    def _send_dac_raw(self):
        addr = self.dac_addr_var.get().strip()
        code = self.dac_raw_code_var.get().strip()
        if not addr:
            messagebox.showwarning("DAC Raw", "Enter/select a DAC I²C address (e.g., 0x4C).")
            return
        if not code.isdigit():
            messagebox.showwarning("DAC Raw", "Raw code must be an integer (0–4095).")
            return
        self._raw_send(f"dac raw {addr} {code}")

    def _raw_send(self, cmdline: str):
        """Send without changing the 'last_command' CSV field (used internally)."""
        if not self.serial_worker:
            return
        try:
            self.serial_worker.send((cmdline.strip() + "\n").encode() if isinstance(cmdline, bytes) else (cmdline.strip() + "\n"))
            self._append_console(f"> {cmdline.strip()}\n")
        except Exception as e:
            self._append_console(f"[ERROR] send failed: {e}\n")

    # --- Parsing + logging ---
    def _handle_line(self, raw: str):
                # ----- Autorun begin/end -----
        m = PATTERNS["autorun_begin"].match(raw)
        if m:
            self.autorun_active = True
            self.autorun_dv.clear()
            self.autorun_dv_a0.clear()
            self.autorun_dv_a1.clear()
            self._a1_dark.clear()
            self._dac_done_for_channel.clear()
            self._last_a0_seen.clear()
            self.csv.log("autorun_begin", None, None, "", raw)
            self._summary_emitted = False
            return

        m = PATTERNS["autorun_end"].match(raw)
        if m:
            # Print both summaries and reset
            if self._summary_emitted:
                return
            self._summary_emitted = True

            self._emit_autorun_summary()
            self.autorun_active = False
            self._dac_done_for_channel.clear()
            self._last_a0_seen.clear()
            self._a1_dark.clear()
            self.csv.log("autorun_end", None, None, "", raw)
            return

        self._append_console(raw + "\n")

        # Autorun begin/end


        # Channel marker
        m = PATTERNS["channel"].match(raw)
        if m:
            self.current_channel = int(m.group(1))
            self.csv.log("channel", self.current_channel, "", "", raw)
            if self.current_channel is not None:
                self._last_a0_seen.pop(self.current_channel, None)
                self._a1_dark.pop(self.current_channel, None)  # clear A1 dark per channel
            return


        # TEMP (on)
        m = PATTERNS["temp_on"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("TEMP_on", self.current_channel, val, "C", raw)
            return

        # TEMP
        m = PATTERNS["temp"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("TEMP", self.current_channel, val, "C", raw)
            return

        # A0 flavors
        # A0 flavors
        for key, event in [("a0_dark","A0_dark"), ("a0_on","A0_on"), ("a0_post","A0_post"), ("a0","A0")]:
            m = PATTERNS[key].match(raw)
            if m:
                val = self._to_float_or_blank(m.group(1))
                self.csv.log(event, self.current_channel, val, "V", raw)
                if self.autorun_active and self.current_channel is not None and isinstance(val, float):
                    self._last_a0_seen[self.current_channel] = val
                    if key == "a0_dark":
                        self._maybe_set_dac_from_a0_dark(self.current_channel, val)
                return
        # A1 flavors
        for key, event in [("a1_dark","A1_dark"), ("a1_on","A1_on"), ("a1_post","A1_post")]:
            m = PATTERNS[key].match(raw)
            if m:
                val = self._to_float_or_blank(m.group(1))
                self.csv.log(event, self.current_channel, val, "V", raw)
                if self.autorun_active and self.current_channel is not None and isinstance(val, float):
                    if key == "a1_dark":
                        self._a1_dark[self.current_channel] = val
                    elif key == "a1_on":
                        dark = self._a1_dark.get(self.current_channel, "")
                        if isinstance(dark, float):
                            self.autorun_dv_a1[self.current_channel] = (val - dark)
                return


        # A1
        m = PATTERNS["a1"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("A1", self.current_channel, val, "V", raw)
            return


                # ΔV (this is A0 / pre-diff from firmware)
        m = PATTERNS["dv"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("deltaV", self.current_channel, val, "V", raw)
            if self.autorun_active and self.current_channel is not None and isinstance(val, float):
                self.autorun_dv[self.current_channel] = val     # keep legacy
                self.autorun_dv_a0[self.current_channel] = val  # new explicit dict
            return


        # I2C scan lines (strict)
        m = PATTERNS["i2c_found"].match(raw)
        if m:
            addr = f"0x{m.group(2).upper()}"
            name_hint = (m.group(1) or "").strip()
            name = I2C_NAMES.get(addr, name_hint or "I2C device")
            self._record_i2c_device(addr, name)
            self.csv.log("i2c_device", None, addr, name, raw)
            return

        # I2C scan lines (flex)
        m = PATTERNS["i2c_found_flex"].search(raw)
        if m:
            addr = f"0x{m.group(1).upper()}"
            name = I2C_NAMES.get(addr, "I2C device")
            self._record_i2c_device(addr, name)
            self.csv.log("i2c_device", None, addr, name, raw)
            return

        # DAC acks
        m = PATTERNS["dac_ack"].match(raw)
        if m:
            v = self._to_float_or_blank(m.group(1))
            self.csv.log("dac_set_volts", None, v, "V", raw)
            return
        m = PATTERNS["dac_raw_ack"].match(raw)
        if m:
            addr = f"0x{m.group(1).upper()}"
            code = m.group(2)
            self.csv.log("dac_set_raw", None, f"{addr}:{code}", "", raw)
            return

        # Warnings/errors
        if raw.startswith("WARNING") or raw.startswith("[ERROR]"):
            self.csv.log("message", self.current_channel, None, "", raw)

    def _record_i2c_device(self, addr: str, name: str):
        self.devices_seen[addr] = name
        self._refresh_devices_panel()
        if name.lower().startswith("dac"):
            dacs = [a for a,n in self.devices_seen.items() if n.lower().startswith("dac")]
            self.detected_dac_var.set(", ".join(sorted(dacs)))

    def _refresh_devices_panel(self):
        for row in self.devs_tree.get_children():
            self.devs_tree.delete(row)
        for addr in sorted(self.devices_seen.keys()):
            self.devs_tree.insert("", END, values=(addr, self.devices_seen[addr]))

    def _emit_autorun_summary(self):
        # A0 (pre-diff) — from firmware ΔV lines
        if not self.autorun_dv_a0:
            self._append_console("[SUMMARY] No ΔV values captured.\n")
            self.csv.log("summary", None, "none", "", "no deltaV captured")
        else:
            self._append_console("\n[SUMMARY] Full Cycle Run — AC Change per LED (A0: pre-diff amp)\n")
            for ch in sorted(self.autorun_dv_a0.keys()):
                dv_v = self.autorun_dv_a0[ch]
                dv_mv_str = self._format_mv(dv_v)
                line = f"LED {ch} AC Change : {dv_mv_str}"
                self._append_console(line + "\n")
                self.csv.log("summary_ac_change_A0", ch, round(dv_v * 1000.0, 3), "mV", line)
            self._append_console("\n")

        # A1 (post-diff) — computed here as A1_on - A1_dark
        if self.autorun_dv_a1:
            self._append_console("[SUMMARY] Full Cycle Run — AC Change per LED (A1: post-diff amp)\n")
            for ch in sorted(self.autorun_dv_a1.keys()):
                dv_v = self.autorun_dv_a1[ch]
                dv_mv_str = self._format_mv(dv_v)
                line = f"LED {ch} AC Change : {dv_mv_str}"
                self._append_console(line + "\n")
                self.csv.log("summary_ac_change_A1", ch, round(dv_v * 1000.0, 3), "mV", line)
            self._append_console("\n")


    # ---- New helper: drive DAC from A0 ----
    def _maybe_set_dac_from_a0_dark(self, channel: int, a0_dark_volts: float):
        """Once per channel during autorun: set DAC = A0_dark - user_ref_margin (clamped)."""
        if channel in self._dac_done_for_channel:
            return
        # Read user margin
        try:
            user_margin = float(self.ref_margin_var.get().strip())
        except Exception:
            user_margin = 0.5  # fallback
        target = a0_dark_volts - user_margin
        if target < DAC_MIN_V: target = DAC_MIN_V
        if target > DAC_MAX_V: target = DAC_MAX_V

        addr = self.dac_addr_var.get().strip()
        cmd = f"dac {addr} {target:.3f}" if addr else f"dac {target:.3f}"
        self._raw_send(cmd)

        self._dac_done_for_channel.add(channel)
        self.csv.log("dac_set_for_channel", channel, round(target, 4), "V",
                     f"DAC set to {target:.3f} V for channel {channel} (A0_dark={a0_dark_volts:.3f} V)")


    # ---- Utils ----
    @staticmethod
    def _to_float_or_blank(s: str):
        try:
            if s.upper() == "N/A" or s.lower() == "nan":
                return ""
            return float(s)
        except Exception:
            return ""

    @staticmethod
    def _format_mv(v_volts: float) -> str:
        try:
            return f"{v_volts*1000:.0f} mV"
        except Exception:
            return "N/A"

def main():
    root = Tk()
    try:
        style = ttk.Style()
        if "vista" in style.theme_names():
            style.theme_use("vista")
        elif "clam" in style.theme_names():
            style.theme_use("clam")
    except Exception:
        pass
    App(root)
    root.mainloop()

if __name__ == "__main__":
    main()
