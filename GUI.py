import os
import sys
import csv
import re
import queue
import threading
from datetime import datetime
from pathlib import Path
from tkinter import Tk, ttk, StringVar, END, DISABLED, N, S, E, W, messagebox
from tkinter import scrolledtext

# ====== User-tweakables ======
BAUD = 115200
READ_TIMEOUT = 0.1  # seconds
GUI_POLL_MS = 50    # ms
SESSIONS_DIRNAME = "PD_Sessions"  # Desktop folder for CSVs

# ====== Known names (helpful labels) ======
I2C_NAMES = {
    "0x40": "TLC59108 LED Driver",
    "0x48": "MAX30205 Temp",
    "0x49": "MAX30205 Temp",
    "0x4A": "MAX30205 Temp",
    "0x4B": "MAX30205 Temp",
    "0x4C": "MAX30205 Temp",
    "0x4D": "MAX30205 Temp",
    "0x4E": "MAX30205 Temp",
    "0x4F": "MAX30205 Temp",
}

# ====== Parsing patterns (robust) ======
PATTERNS = {
    "autorun_begin": re.compile(r"^===\s*(?:Full\s*Cycle\s*Run|Auto\-tuned\s*Autorun)", re.I),
    "autorun_end":   re.compile(r"^\[INFO\]\s*Auto\-tuned\s*autorun\s*complete\.|^===\s*Autorun\s*complete\s*===\s*$", re.I),
    "channel":       re.compile(r"^===\s*Channel\s+(\d+)\b.*$", re.I),

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

    "gain":          re.compile(r"^Gain\s*\(\s*ΔA1/ΔA0\s*\)\s*:\s*([\-NaN\d\.]+)", re.I),
    "dv":            re.compile(r"^ΔV\s*=\s*([\-NaN\d\.]+)", re.I),

    "i2c_found":      re.compile(r"^(.*)\s:\s0x([0-9A-Fa-f]{2})$"),
    "i2c_found_flex": re.compile(r"(?:found.*at\s*)?0x([0-9A-Fa-f]{2})", re.I),

    # gslope-specific
    "gslope_begin":    re.compile(r"^===\s*gslope\s*ch\s+(\d+)\s*===", re.I),
    "gslope_point":    re.compile(r"^PWM\s+(\d+).*?A0_on:\s*([\-.\d]+).*?A1_on:\s*([\-.\d]+).*?dA0:\s*([\-.\d]+).*?dA1:\s*([\-.\d]+)", re.I),
    "gslope_result":   re.compile(r"^gslope\s*result\s*—\s*slope:\s*([\-.\d]+)\s*,\s*intercept:\s*([\-.\d]+)\s*,\s*R\^2:\s*([\-.\d]+)", re.I),
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
            if isinstance(s, str):
                self.ser.write(s.encode())
            else:
                self.ser.write(s)

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
        self.root.title("Photodiode Board Controller — Serial GUI (No DAC)")
        self.root.geometry("1120x740")

        # State
        self.serial_worker = None
        self.line_queue = queue.Queue()
        self.csv = CsvLogger()
        self.current_channel = None

        # Autorun aggregation state
        self.autorun_active = False
        self._summary_emitted = False
        self.autorun_dv_a0 = {}     # A0 ΔV
        self.autorun_dv_a1 = {}     # A1 ΔV
        self.autorun_gain  = {}     # ΔA1/ΔA0
        self._a1_dark = {}          # channel -> last A1_dark

        # gslope state (for parsing)
        self._gslope_current_ch = None

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

        self._mk_button(btns, "Full Cycle Run (r)",   lambda: self.send_cmd("r"),   0, 0)
        self._mk_button(btns, "Auto-tuned ΔA0 (r0t)", lambda: self.send_cmd("r0t"), 0, 1)
        self._mk_button(btns, "Temp (t)",             lambda: self.send_cmd("t"),   0, 2)
        self._mk_button(btns, "A0 (a0)",              lambda: self.send_cmd("a0"),  0, 3)
        self._mk_button(btns, "A1 (a1)",              lambda: self.send_cmd("a1"),  0, 4)
        self._mk_button(btns, "I²C Scan",             self.scan_i2c,                0, 5)
        self._mk_button(btns, "Help/menu",            lambda: self.send_cmd("help"),0, 6)

        # ===== gslope controls =====
        gs = ttk.LabelFrame(btns, text="gslope (linear fit ΔA1 vs ΔA0)")
        gs.grid(row=3, column=0, columnspan=7, sticky=E+W, pady=(10,0))

        ttk.Label(gs, text="Channel:").grid(row=0, column=0, padx=(4,4), sticky=E)
        self.gs_ch_var = StringVar(value="all")
        self.gs_ch_combo = ttk.Combobox(gs, textvariable=self.gs_ch_var, width=6, state="readonly")
        self.gs_ch_combo["values"] = ["all","0","1","2","3"]
        self.gs_ch_combo.grid(row=0, column=1, padx=(0,12), sticky=W)

        ttk.Label(gs, text="PWM list (CSV, optional):").grid(row=0, column=2, padx=(4,4), sticky=E)
        self.gs_pwm_var = StringVar(value="120,160,200,240")
        self.gs_pwm_entry = ttk.Entry(gs, textvariable=self.gs_pwm_var, width=24)
        self.gs_pwm_entry.grid(row=0, column=3, padx=(0,12), sticky=W)

        self._mk_button(gs, "Run gslope", self._run_gslope_click, 0, 4)

        for c in range(5):
            gs.grid_columnconfigure(c, weight=1)

        # ===== LED Intensity: Textboxes instead of sliders =====
        intensity = ttk.LabelFrame(btns, text="LED Intensity (0–255): press Enter to apply")
        intensity.grid(row=2, column=0, columnspan=7, sticky=E+W, pady=(8,0))

        self._led_vars = []
        for ch in range(4):
            ttk.Label(intensity, text=f"LED {ch}").grid(row=0, column=ch*3, padx=(4,2), sticky=E)
            var = StringVar(value="0")
            ent = ttk.Entry(intensity, textvariable=var, width=6, justify="right")
            ent.grid(row=0, column=ch*3+1, padx=(0,2), sticky=W)
            ent.bind("<Return>", lambda e, c=ch: self._apply_pwm_from_entry(c))
            btn = ttk.Button(intensity, text="Set", command=lambda c=ch: self._apply_pwm_from_entry(c))
            btn.grid(row=0, column=ch*3+2, padx=(0,6), sticky=W)
            self._led_vars.append(var)
        for col in range(12):
            intensity.grid_columnconfigure(col, weight=1)

        # Quick OFF/ON buttons (ON uses 255)
        led_on = ttk.LabelFrame(btns, text="LED ON @255")
        led_on.grid(row=1, column=0, columnspan=3, sticky=E+W, pady=(8,0))
        for ch in range(4):
            self._mk_button(led_on, f"LED {ch} ON",  lambda c=ch: self._set_and_send_pwm(c, 255), 0, ch)

        led_off = ttk.LabelFrame(btns, text="LED OFF")
        led_off.grid(row=1, column=3, columnspan=3, sticky=E+W, pady=(8,0))
        for ch in range(4):
            self._mk_button(led_off, f"LED {ch} OFF", lambda c=ch: self._set_and_send_pwm(c, 0), 0, ch)

        # ===== Devices panel =====
        devs = ttk.LabelFrame(root, text="I²C Devices (from last scan)", padding=8)
        devs.grid(row=2, column=0, sticky=E+W, padx=8, pady=(0,8))
        cols = ("addr", "name")
        self.devs_tree = ttk.Treeview(devs, columns=cols, show="headings", height=6)
        self.devs_tree.heading("addr", text="Address")
        self.devs_tree.heading("name", text="Name")
        self.devs_tree.column("addr", width=90, anchor="center")
        self.devs_tree.column("name", width=300, anchor="w")
        self.devs_tree.grid(row=0, column=0, sticky=E+W)
        devs.grid_columnconfigure(0, weight=1)

        # ===== Console =====
        self.console = scrolledtext.ScrolledText(root, wrap="word", height=22)
        self.console.grid(row=3, column=0, sticky=N+S+E+W, padx=8, pady=(0,8))
        self._append_console(self._append_notice)

        # Layout stretch
        root.grid_rowconfigure(3, weight=1)
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
        if not self.serial_worker:
            messagebox.showwarning("Not Connected", "Connect to a serial port first.")
            return

        lc = cmd.strip().lower()
        # Reset autorun aggregation at the start of r or r0t
        if lc in ("r", "r0t"):
            self.autorun_active = True
            self.autorun_dv_a0.clear()
            self.autorun_dv_a1.clear()
            self.autorun_gain.clear()
            self._a1_dark.clear()
            self._summary_emitted = False

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
        self.send_cmd("scan")

    def _run_gslope_click(self):
        """
        Build and send the appropriate 'gslope' command from controls.
        """
        ch = self.gs_ch_var.get().strip().lower()
        pwm_csv = self.gs_pwm_var.get().strip().replace(" ", "")
        # Determine command shape
        if ch == "all" and (pwm_csv == "" or pwm_csv.lower() == "default"):
            cmd = "gslope"
        elif ch == "all" and pwm_csv:
            cmd = f"gslope {pwm_csv}"
        elif ch in ("0","1","2","3") and (pwm_csv == "" or pwm_csv.lower() == "default"):
            cmd = f"gslope {ch}"
        else:
            cmd = f"gslope {ch} {pwm_csv}"
        self.send_cmd(cmd)

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

    # ---- PWM helpers ----
    def _apply_pwm_from_entry(self, ch: int):
        txt = self._led_vars[ch].get().strip()
        if txt == "":
            messagebox.showwarning("PWM", f"Enter a value for LED {ch} (0–255).")
            return
        try:
            val = int(float(txt))
        except ValueError:
            messagebox.showwarning("PWM", f"Invalid value '{txt}' for LED {ch}. Use 0–255.")
            return
        if val < 0 or val > 255:
            messagebox.showwarning("PWM", f"Out of range for LED {ch}: {val}. Use 0–255.")
            return
        self._led_vars[ch].set(str(val))
        self._raw_send(f"pwm {ch} {val}")

    def _set_and_send_pwm(self, ch: int, val: int):
        self._led_vars[ch].set(str(val))
        self._raw_send(f"pwm {ch} {val}")

    # ---- Raw sender ----
    def _raw_send(self, cmdline: str):
        if not self.serial_worker:
            return
        try:
            payload = (cmdline.strip() + "\n")
            self.serial_worker.send(payload)
            self._append_console(f"> {cmdline.strip()}\n")
            self.csv.set_last_command(cmdline.strip())
        except Exception as e:
            self._append_console(f"[ERROR] send failed: {e}\n")

    # ---- Parsing + logging ----
    def _handle_line(self, raw: str):
        # gslope begin?
        m = PATTERNS["gslope_begin"].match(raw)
        if m:
            self._gslope_current_ch = int(m.group(1))
            self.csv.log("gslope_begin", self._gslope_current_ch, None, "", raw)
            self._append_console(raw + "\n")
            return

        # gslope point
        m = PATTERNS["gslope_point"].match(raw)
        if m:
            pwm = m.group(1); a0_on = m.group(2); a1_on = m.group(3); dA0 = m.group(4); dA1 = m.group(5)
            ch = self._gslope_current_ch if self._gslope_current_ch is not None else ""
            self.csv.log("gslope_point_pwm", ch, int(pwm), "", raw)
            self.csv.log("gslope_point_dA0", ch, float(dA0), "V", raw)
            self.csv.log("gslope_point_dA1", ch, float(dA1), "V", raw)
            self._append_console(raw + "\n")
            return

        # gslope result
        m = PATTERNS["gslope_result"].match(raw)
        if m:
            slope = float(m.group(1))
            intercept = float(m.group(2))
            r2 = float(m.group(3))
            ch = self._gslope_current_ch if self._gslope_current_ch is not None else ""
            self.csv.log("gslope_slope", ch, slope, "", raw)
            self.csv.log("gslope_intercept", ch, intercept, "V", raw)
            self.csv.log("gslope_r2", ch, r2, "", raw)
            self._append_console(raw + "\n")
            # do not clear _gslope_current_ch yet; another channel may follow with new header
            return

        # Autorun begin?
        if PATTERNS["autorun_begin"].match(raw):
            self.autorun_active = True
            self.autorun_dv_a0.clear()
            self.autorun_dv_a1.clear()
            self.autorun_gain.clear()
            self._a1_dark.clear()
            self.csv.log("autorun_begin", None, None, "", raw)
            self._summary_emitted = False
            self._append_console(raw + "\n")
            return

        # Autorun end?
        if PATTERNS["autorun_end"].match(raw):
            if not self._summary_emitted:
                self._emit_autorun_summary()
                self._summary_emitted = True
            self.autorun_active = False
            self.csv.log("autorun_end", None, None, "", raw)
            self._append_console(raw + "\n")
            return

        self._append_console(raw + "\n")

        # Channel marker
        m = PATTERNS["channel"].match(raw)
        if m:
            self.current_channel = int(m.group(1))
            self.csv.log("channel", self.current_channel, "", "", raw)
            if self.current_channel is not None:
                self._a1_dark.pop(self.current_channel, None)
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
        for key, event in [("a0_dark","A0_dark"), ("a0_on","A0_on"), ("a0_post","A0_post"), ("a0","A0")]:
            m = PATTERNS[key].match(raw)
            if m:
                val = self._to_float_or_blank(m.group(1))
                self.csv.log(event, self.current_channel, val, "V", raw)
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

        # ΔV (A0 / pre-diff from firmware)
        m = PATTERNS["dv"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("deltaV_A0", self.current_channel, val, "V", raw)
            if self.autorun_active and self.current_channel is not None and isinstance(val, float):
                self.autorun_dv_a0[self.current_channel] = val
            return

        # GAIN line (ΔA1/ΔA0)
        m = PATTERNS["gain"].match(raw)
        if m:
            val = self._to_float_or_blank(m.group(1))
            self.csv.log("gain_delta", self.current_channel, val, "", raw)
            if self.autorun_active and self.current_channel is not None and isinstance(val, float):
                self.autorun_gain[self.current_channel] = val
            return

        # I2C (strict)
        m = PATTERNS["i2c_found"].match(raw)
        if m:
            addr = f"0x{m.group(2).upper()}"
            name_hint = (m.group(1) or "").strip()
            name = I2C_NAMES.get(addr, name_hint or "I2C device")
            self._record_i2c_device(addr, name)
            self.csv.log("i2c_device", None, addr, name, raw)
            return

        # I2C (flex)
        m = PATTERNS["i2c_found_flex"].search(raw)
        if m:
            addr = f"0x{m.group(1).upper()}"
            name = I2C_NAMES.get(addr, "I2C device")
            self._record_i2c_device(addr, name)
            self.csv.log("i2c_device", None, addr, name, raw)
            return

        # Warnings/errors
        if raw.startswith("WARNING") or raw.startswith("[ERROR]"):
            self.csv.log("message", self.current_channel, None, "", raw)

    def _record_i2c_device(self, addr: str, name: str):
        self.devices_seen[addr] = name
        self._refresh_devices_panel()

    def _refresh_devices_panel(self):
        for row in self.devs_tree.get_children():
            self.devs_tree.delete(row)
        for addr in sorted(self.devices_seen.keys()):
            self.devs_tree.insert("", END, values=(addr, self.devices_seen[addr]))

    def _emit_autorun_summary(self):
        # A0 summary
        if not self.autorun_dv_a0:
            self._append_console("[SUMMARY] No ΔV(A0) values captured.\n")
            self.csv.log("summary", None, "none", "", "no deltaV captured")
        else:
            self._append_console("\n[SUMMARY] Full Cycle Run — AC Change per LED (A0: pre-diff amp)\n")
            for ch in sorted(self.autorun_dv_a0.keys()):
                dv_v = self.autorun_dv_a0[ch]
                dv_mv = f"{dv_v*1000:.0f} mV"
                line = f"LED {ch} AC Change : {dv_mv}"
                self._append_console(line + "\n")
                self.csv.log("summary_ac_change_A0", ch, round(dv_v * 1000.0, 3), "mV", line)
            self._append_console("\n")

        # A1 summary
        if self.autorun_dv_a1:
            self._append_console("[SUMMARY] Full Cycle Run — AC Change per LED (A1: post-diff amp)\n")
            for ch in sorted(self.autorun_dv_a1.keys()):
                dv_v = self.autorun_dv_a1[ch]
                dv_mv = f"{dv_v*1000:.0f} mV"
                line = f"LED {ch} AC Change : {dv_mv}"
                self._append_console(line + "\n")
                self.csv.log("summary_ac_change_A1", ch, round(dv_v * 1000.0, 3), "mV", line)
            self._append_console("\n")

        # Gain summary (ΔA1/ΔA0)
        if self.autorun_gain:
            self._append_console("[SUMMARY] Full Cycle Run — Gain per LED (ΔA1/ΔA0)\n")
            for ch in sorted(self.autorun_gain.keys()):
                g = self.autorun_gain[ch]
                line = f"LED {ch} Gain : {g:.3f}"
                self._append_console(line + "\n")
                self.csv.log("summary_gain", ch, g, "", line)
            self._append_console("\n")

    # ---- Utils ----
    @staticmethod
    def _to_float_or_blank(s: str):
        try:
            if s.upper() == "N/A" or s.lower() == "nan":
                return ""
            return float(s)
        except Exception:
            return ""

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
