import os
import re
import csv
import time
import math
import queue
import threading
from datetime import datetime

import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.scrolledtext import ScrolledText

import serial
import serial.tools.list_ports

# ============================
# USER CONFIG: channel mapping
# ============================
# YOUR confirmed mapping:
#   CH1 = 1450 close
#   CH2 = 1650 far
#   CH3 = 1450 far
#   CH4 = 1650 close
CH_1450_CLOSE = 1
CH_1450_FAR   = 3
CH_1650_CLOSE = 4
CH_1650_FAR   = 2

# ============================
# Firmware line matching
# ============================
# RX line example:
# RX CH2: I=45.00 mA, A0_avg=1234.56 counts, 1.234567 V
RX_RE = re.compile(r"^RX\s+CH(\d+):.*?([0-9]+(?:\.[0-9]+)?)\s*V\s*$")

# A0 line example:
# A0_avg(50) = 1234.56 counts, 1.234567 V
A0_RE = re.compile(r"^A0_avg\(\d+\)\s*=\s*.*?,\s*([0-9]+(?:\.[0-9]+)?)\s*V\s*$")


def desktop_log_folder() -> str:
    home = os.path.expanduser("~")
    desktop = os.path.join(home, "Desktop")
    base = desktop if os.path.isdir(desktop) else home
    folder = os.path.join(base, "MPQ3326_Logs")
    os.makedirs(folder, exist_ok=True)
    return folder


def now_iso() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S.%f%z")


def sanitize_filename(name: str) -> str:
    # Remove characters illegal on Windows filenames: \ / : * ? " < > |
    illegal = r'\/:*?"<>|'
    out = "".join("_" if c in illegal else c for c in name.strip())
    out = out.replace("\n", "_").replace("\r", "_")
    return out if out else "run"


def safe_log(x: float):
    if x is None or x <= 0:
        return ""
    return math.log(x)


def safe_ratio(a: float, b: float):
    if a is None or b is None:
        return ""
    if b == 0:
        return ""
    return a / b


def safe_log_ratio(a: float, b: float):
    if a is None or b is None:
        return ""
    if a <= 0 or b <= 0:
        return ""
    return math.log(a / b)


class SerialReader(threading.Thread):
    def __init__(self, ser: serial.Serial, out_q: queue.Queue, stop_evt: threading.Event):
        super().__init__(daemon=True)
        self.ser = ser
        self.out_q = out_q
        self.stop_evt = stop_evt

    def run(self):
        buf = b""
        while not self.stop_evt.is_set():
            try:
                chunk = self.ser.read(256)
                if chunk:
                    buf += chunk
                    while b"\n" in buf:
                        line, buf = buf.split(b"\n", 1)
                        line = line.replace(b"\r", b"").decode(errors="replace").strip()
                        if line:
                            self.out_q.put(line)
                else:
                    time.sleep(0.01)
            except Exception as e:
                self.out_q.put(f"[SERIAL_ERROR] {e}")
                break


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("MPQ3326 Serial GUI")

        # Serial
        self.ser = None
        self.reader_stop = threading.Event()
        self.reader_thread = None
        self.rx_q = queue.Queue()

        # Continuous mode state machine
        self.cont_running = False
        self.stage = "idle"  # idle, awaiting_dark, awaiting_r1..awaiting_r4
        self.dark_v = None
        self.cycle_index = 0
        self.row = {}

        # CSV
        self.csv_file = None
        self.csv_writer = None
        self.csv_path = None

        # UI vars
        self.port_var = tk.StringVar()
        self.baud_var = tk.StringVar(value="115200")
        self.dark_samples_var = tk.StringVar(value="50")   # a0 N (default 50)
        self.cycle_delay_var = tk.StringVar(value="0.20")  # seconds between cycles
        self.log_to_csv_var = tk.BooleanVar(value=True)
        self.run_label_var = tk.StringVar(value="70%_3days")  # user edits manually

        self._build_ui()
        self._refresh_ports()
        self.root.after(50, self._poll_serial_queue)

    def _build_ui(self):
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="x")

        ttk.Label(top, text="Port:").pack(side="left")
        self.port_combo = ttk.Combobox(top, textvariable=self.port_var, width=18, state="readonly")
        self.port_combo.pack(side="left", padx=(4, 10))
        ttk.Button(top, text="Refresh", command=self._refresh_ports).pack(side="left")

        ttk.Label(top, text="Baud:").pack(side="left", padx=(10, 4))
        ttk.Entry(top, textvariable=self.baud_var, width=8).pack(side="left")

        self.btn_connect = ttk.Button(top, text="Connect", command=self.connect)
        self.btn_connect.pack(side="left", padx=(10, 4))
        self.btn_disconnect = ttk.Button(top, text="Disconnect", command=self.disconnect, state="disabled")
        self.btn_disconnect.pack(side="left")

        mid = ttk.Frame(self.root, padding=8)
        mid.pack(fill="x")

        ttk.Label(mid, text="Command:").pack(side="left")
        self.cmd_entry = ttk.Entry(mid, width=40)
        self.cmd_entry.pack(side="left", padx=(6, 6))
        self.cmd_entry.bind("<Return>", lambda e: self.send_command())
        ttk.Button(mid, text="Send", command=self.send_command).pack(side="left")

        quick = ttk.Frame(self.root, padding=(8, 0, 8, 8))
        quick.pack(fill="x")

        ttk.Button(quick, text="help", command=lambda: self._send("help")).pack(side="left")
        ttk.Button(quick, text="dump", command=lambda: self._send("dump")).pack(side="left", padx=4)
        ttk.Button(quick, text="off", command=lambda: self._send("off")).pack(side="left")

        ttk.Separator(quick, orient="vertical").pack(side="left", fill="y", padx=10)

        ttk.Button(quick, text="r1", command=lambda: self._send("r1")).pack(side="left")
        ttk.Button(quick, text="r2", command=lambda: self._send("r2")).pack(side="left", padx=4)
        ttk.Button(quick, text="r3", command=lambda: self._send("r3")).pack(side="left")
        ttk.Button(quick, text="r4", command=lambda: self._send("r4")).pack(side="left", padx=4)

        ttk.Separator(quick, orient="vertical").pack(side="left", fill="y", padx=10)

        ttk.Button(quick, text="a0 50", command=lambda: self._send("a0 50")).pack(side="left")

        box = ttk.LabelFrame(self.root, text="Continuous Mode (cycles r1→r2→r3→r4 until Stop)", padding=8)
        box.pack(fill="x", padx=8, pady=(0, 8))

        ttk.Label(box, text="Run label (used in filename + CSV):").grid(row=0, column=0, sticky="w")
        ttk.Entry(box, textvariable=self.run_label_var, width=20).grid(row=0, column=1, sticky="w", padx=(6, 12))

        ttk.Label(box, text="Dark samples (a0 N):").grid(row=0, column=2, sticky="w")
        ttk.Entry(box, textvariable=self.dark_samples_var, width=6).grid(row=0, column=3, sticky="w", padx=(6, 12))

        ttk.Label(box, text="Delay between cycles (s):").grid(row=0, column=4, sticky="w")
        ttk.Entry(box, textvariable=self.cycle_delay_var, width=6).grid(row=0, column=5, sticky="w", padx=(6, 12))

        ttk.Checkbutton(box, text="Log to CSV on Desktop", variable=self.log_to_csv_var).grid(
            row=0, column=6, sticky="w", padx=(6, 6)
        )

        self.btn_start = ttk.Button(box, text="Start Continuous", command=self.start_continuous, state="disabled")
        self.btn_start.grid(row=0, column=7, sticky="w", padx=(6, 6))

        self.btn_stop = ttk.Button(box, text="Stop", command=self.stop_continuous, state="disabled")
        self.btn_stop.grid(row=0, column=8, sticky="w")

        self.status = ttk.Label(box, text="Idle")
        self.status.grid(row=1, column=0, columnspan=9, sticky="w", pady=(6, 0))

        logframe = ttk.Frame(self.root, padding=8)
        logframe.pack(fill="both", expand=True)

        self.log = ScrolledText(logframe, height=18, wrap="word")
        self.log.pack(fill="both", expand=True)

    def _refresh_ports(self):
        ports = [p.device for p in serial.tools.list_ports.comports()]
        self.port_combo["values"] = ports
        if ports and not self.port_var.get():
            self.port_var.set(ports[0])

    def _log(self, msg: str):
        self.log.insert("end", msg + "\n")
        self.log.see("end")

    def connect(self):
        if self.ser:
            return
        port = self.port_var.get().strip()
        if not port:
            messagebox.showerror("Error", "Select a COM port.")
            return
        try:
            baud = int(self.baud_var.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Invalid baud rate.")
            return

        try:
            self.ser = serial.Serial(port, baudrate=baud, timeout=0.05)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open {port}: {e}")
            return

        self.reader_stop.clear()
        self.reader_thread = SerialReader(self.ser, self.rx_q, self.reader_stop)
        self.reader_thread.start()

        self._log(f"[INFO] Connected to {port} @ {baud}")
        self.btn_connect.configure(state="disabled")
        self.btn_disconnect.configure(state="normal")
        self.btn_start.configure(state="normal")

    def disconnect(self):
        self.stop_continuous()
        if not self.ser:
            return
        try:
            self.reader_stop.set()
            time.sleep(0.05)
            self.ser.close()
        except Exception:
            pass
        self.ser = None
        self._log("[INFO] Disconnected.")
        self.btn_connect.configure(state="normal")
        self.btn_disconnect.configure(state="disabled")
        self.btn_start.configure(state="disabled")

    def send_command(self):
        cmd = self.cmd_entry.get().strip()
        if not cmd:
            return
        self._send(cmd)
        self.cmd_entry.delete(0, "end")

    def _send(self, cmd: str):
        if not self.ser:
            self._log("[WARN] Not connected.")
            return
        try:
            self.ser.write((cmd + "\n").encode())
            self._log(f"> {cmd}")
        except Exception as e:
            self._log(f"[SERIAL_WRITE_ERROR] {e}")

    # ----------------------------
    # Continuous Mode
    # ----------------------------
    def start_continuous(self):
        if not self.ser:
            messagebox.showerror("Error", "Connect to serial first.")
            return
        if self.cont_running:
            return

        # Parse/validate dark samples
        try:
            dark_n = int(self.dark_samples_var.get().strip())
            if dark_n <= 0:
                raise ValueError()
        except ValueError:
            messagebox.showerror("Error", "Dark samples must be a positive integer.")
            return

        run_label = self.run_label_var.get().strip()
        if not run_label:
            messagebox.showerror("Error", "Run label is required (e.g., 70%_3days).")
            return

        # Open CSV if enabled
        if self.log_to_csv_var.get():
            folder = desktop_log_folder()
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            prefix = sanitize_filename(run_label)
            self.csv_path = os.path.join(folder, f"{prefix}_{ts}.csv")

            try:
                self.csv_file = open(self.csv_path, "w", newline="", encoding="utf-8")
                self.csv_writer = csv.writer(self.csv_file)

                # CSV header (raw + derived features)
                self.csv_writer.writerow([
                    "index",
                    "timestamp_iso",
                    "run_label",
                    "dark_V",

                    # Raw voltages by channel
                    "CH1_V", "CH2_V", "CH3_V", "CH4_V",

                    # Dark-subtracted
                    "CH1_ds", "CH2_ds", "CH3_ds", "CH4_ds",

                    # Log-space of dark-subtracted (blank if <=0)
                    "ln_CH1_ds", "ln_CH2_ds", "ln_CH3_ds", "ln_CH4_ds",

                    # Meaningful ratios on dark-subtracted (based on your mapping)
                    "ratio_1450_far_over_close",        # (CH3_ds / CH1_ds)
                    "ratio_1650_far_over_close",        # (CH2_ds / CH4_ds)
                    "ratio_1450_close_over_1650_close", # (CH1_ds / CH4_ds)
                    "ratio_1450_far_over_1650_far",     # (CH3_ds / CH2_ds)

                    # Log ratios (blank if operands <=0)
                    "ln_ratio_1450_far_over_close",
                    "ln_ratio_1650_far_over_close",
                    "ln_ratio_1450_close_over_1650_close",
                    "ln_ratio_1450_far_over_1650_far",
                ])
                self.csv_file.flush()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to open CSV:\n{e}")
                self._close_csv()
                return

            self._log(f"[CONT] CSV logging enabled: {self.csv_path}")
        else:
            self._log("[CONT] CSV logging disabled.")

        self.cont_running = True
        self.stage = "awaiting_dark"
        self.dark_v = None
        self.cycle_index = 0
        self.row = {}

        self.btn_start.configure(state="disabled")
        self.btn_stop.configure(state="normal")

        self.status.configure(text="Starting... Taking dark (a0 N) then cycling r1→r4 continuously.")
        self._send("off")
        self.root.after(150, lambda: self._send(f"a0 {dark_n}"))

    def stop_continuous(self):
        if not self.cont_running:
            return
        self.cont_running = False
        self.stage = "idle"
        self.status.configure(text="Stopping...")
        if self.ser:
            self._send("off")
        self._close_csv()
        self.btn_start.configure(state="normal" if self.ser else "disabled")
        self.btn_stop.configure(state="disabled")
        self.status.configure(text="Idle")
        self._log("[CONT] Stopped.")

    def _close_csv(self):
        try:
            if self.csv_file:
                self.csv_file.flush()
                self.csv_file.close()
        except Exception:
            pass
        self.csv_file = None
        self.csv_writer = None

    def _start_cycle(self):
        self.cycle_index += 1
        self.row = {
            "timestamp_iso": now_iso(),
            "CH1_V": None,
            "CH2_V": None,
            "CH3_V": None,
            "CH4_V": None,
        }
        self.stage = "awaiting_r1"
        self.status.configure(text=f"Cycle {self.cycle_index}: running r1→r2→r3→r4")
        self._send("r1")

    def _compute_and_write_row(self):
        run_label = self.run_label_var.get().strip()

        # Raw voltages
        ch = {
            1: self.row["CH1_V"],
            2: self.row["CH2_V"],
            3: self.row["CH3_V"],
            4: self.row["CH4_V"],
        }

        # Dark subtraction
        ds = {}
        for k in [1, 2, 3, 4]:
            if ch[k] is None or self.dark_v is None:
                ds[k] = None
            else:
                ds[k] = ch[k] - self.dark_v

        # Logs of dark-subtracted
        ln_ds = {k: safe_log(ds[k]) if ds[k] is not None else "" for k in [1, 2, 3, 4]}

        # Mapped signals for ratios
        v1450c = ds.get(CH_1450_CLOSE)  # CH1_ds
        v1450f = ds.get(CH_1450_FAR)    # CH3_ds
        v1650c = ds.get(CH_1650_CLOSE)  # CH4_ds
        v1650f = ds.get(CH_1650_FAR)    # CH2_ds

        # Ratios (dark-subtracted domain)
        ratio_1450_far_close = safe_ratio(v1450f, v1450c)  # CH3/CH1
        ratio_1650_far_close = safe_ratio(v1650f, v1650c)  # CH2/CH4
        ratio_1450c_1650c    = safe_ratio(v1450c, v1650c)  # CH1/CH4
        ratio_1450f_1650f    = safe_ratio(v1450f, v1650f)  # CH3/CH2

        # Log ratios
        ln_ratio_1450_far_close = safe_log_ratio(v1450f, v1450c)
        ln_ratio_1650_far_close = safe_log_ratio(v1650f, v1650c)
        ln_ratio_1450c_1650c    = safe_log_ratio(v1450c, v1650c)
        ln_ratio_1450f_1650f    = safe_log_ratio(v1450f, v1650f)

        def fmt(x):
            if x is None or x == "":
                return ""
            if isinstance(x, float):
                return f"{x:.12g}"
            return str(x)

        if self.csv_writer and self.csv_file:
            self.csv_writer.writerow([
                self.cycle_index,
                self.row["timestamp_iso"],
                run_label,
                fmt(self.dark_v),

                fmt(ch[1]), fmt(ch[2]), fmt(ch[3]), fmt(ch[4]),
                fmt(ds[1]), fmt(ds[2]), fmt(ds[3]), fmt(ds[4]),
                fmt(ln_ds[1]), fmt(ln_ds[2]), fmt(ln_ds[3]), fmt(ln_ds[4]),

                fmt(ratio_1450_far_close),
                fmt(ratio_1650_far_close),
                fmt(ratio_1450c_1650c),
                fmt(ratio_1450f_1650f),

                fmt(ln_ratio_1450_far_close),
                fmt(ln_ratio_1650_far_close),
                fmt(ln_ratio_1450c_1650c),
                fmt(ln_ratio_1450f_1650f),
            ])
            self.csv_file.flush()

    def _finish_cycle(self):
        self._compute_and_write_row()

        # Status (guard against None)
        def s(v):
            return "NA" if v is None else f"{v:.6f}"

        self.status.configure(
            text=f"Cycle {self.cycle_index} logged | Dark={s(self.dark_v)} V | "
                 f"CH1={s(self.row['CH1_V'])} CH2={s(self.row['CH2_V'])} "
                 f"CH3={s(self.row['CH3_V'])} CH4={s(self.row['CH4_V'])}"
        )

        # Next cycle
        try:
            delay_s = float(self.cycle_delay_var.get().strip())
            if delay_s < 0:
                delay_s = 0.0
        except ValueError:
            delay_s = 0.2

        if self.cont_running:
            self.root.after(int(delay_s * 1000), self._start_cycle)

    # ----------------------------
    # Serial queue processing
    # ----------------------------
    def _poll_serial_queue(self):
        try:
            while True:
                line = self.rx_q.get_nowait()
                self._handle_serial_line(line)
        except queue.Empty:
            pass
        self.root.after(50, self._poll_serial_queue)

    def _handle_serial_line(self, line: str):
        self._log(line)

        if not self.cont_running:
            return

        # Dark captured once at start
        if self.stage == "awaiting_dark":
            m = A0_RE.match(line)
            if m:
                self.dark_v = float(m.group(1))
                self._log(f"[CONT] Dark captured: {self.dark_v:.6f} V")
                self._start_cycle()
            return

        m = RX_RE.match(line)
        if not m:
            return

        ch = int(m.group(1))
        v = float(m.group(2))

        if self.stage == "awaiting_r1" and ch == 1:
            self.row["CH1_V"] = v
            self.stage = "awaiting_r2"
            self._send("r2")
            return

        if self.stage == "awaiting_r2" and ch == 2:
            self.row["CH2_V"] = v
            self.stage = "awaiting_r3"
            self._send("r3")
            return

        if self.stage == "awaiting_r3" and ch == 3:
            self.row["CH3_V"] = v
            self.stage = "awaiting_r4"
            self._send("r4")
            return

        if self.stage == "awaiting_r4" and ch == 4:
            self.row["CH4_V"] = v
            self._finish_cycle()
            return


def main():
    root = tk.Tk()
    app = App(root)

    def on_close():
        try:
            app.disconnect()
        except Exception:
            pass
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.mainloop()


if __name__ == "__main__":
    main()
