import os
import sys
import csv
import re
import queue
import threading
import time
from datetime import datetime
from pathlib import Path
from tkinter import Tk, ttk, StringVar, END, DISABLED, NORMAL, N, S, E, W, messagebox
from tkinter import scrolledtext
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ====== User-tweakables ======
BAUD = 115200
READ_TIMEOUT = 0.1
GUI_POLL_MS = 20    # Faster polling for smooth ramp
SESSIONS_DIRNAME = "PD_Sessions"

# ====== Parsing Patterns ======
PATTERNS = {
    "temp": re.compile(r"TEMP.*:\s*([-\d\.]+)"),
    "a0":   re.compile(r"A0.*:\s*([-\d\.]+)"),
    "a1":   re.compile(r"A1.*:\s*([-\d\.]+)"),
    "a0_dark": re.compile(r"A0 \(dark\):\s*([-\d\.]+)"),
    "a1_dark": re.compile(r"A1 \(dark\):\s*([-\d\.]+)"),
    "a0_on":   re.compile(r"A0 \(on\):\s*([-\d\.]+)"),
    "a1_on":   re.compile(r"A1 \(on\):\s*([-\d\.]+)")
}

class SerialReader(threading.Thread):
    def __init__(self, serial_port, data_queue):
        super().__init__()
        self.serial = serial_port
        self.queue = data_queue
        self.running = True
        self.daemon = True

    def run(self):
        while self.running:
            try:
                if self.serial.in_waiting:
                    line = self.serial.readline()
                    if line:
                        try:
                            text = line.decode('utf-8', errors='ignore').strip()
                            self.queue.put(text)
                        except:
                            pass
                else:
                    time.sleep(0.01)
            except Exception:
                self.running = False

    def stop(self):
        self.running = False

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("PD Board Control v2 (Ramp & Fast)")
        
        # Data & State
        self.ser = None
        self.reader = None
        self.msg_queue = queue.Queue()
        self.is_connected = False
        
        # Ramp State
        self.ramp_running = False
        self.ramp_data = [] # List of (pwm, a0, a1)
        self.ramp_channel = 0
        self.ramp_step = 10
        self.ramp_delay = 0.05 # Delay between step and measure
        
        # UI Setup
        self._setup_ui()
        
        # Start GUI Loop
        self.root.after(GUI_POLL_MS, self._gui_loop)

    def _setup_ui(self):
        # Layout: Left Control Panel, Right Log/Plot
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill="both", expand=True)
        
        # --- Connection Bar ---
        conn_frame = ttk.LabelFrame(main_frame, text="Connection")
        conn_frame.pack(fill="x", pady=5)
        
        ttk.Label(conn_frame, text="Port:").pack(side="left", padx=5)
        self.port_var = StringVar(value="COM4")
        ttk.Entry(conn_frame, textvariable=self.port_var, width=10).pack(side="left")
        
        self.btn_connect = ttk.Button(conn_frame, text="Connect", command=self._toggle_connect)
        self.btn_connect.pack(side="left", padx=10)
        
        # --- Tabs ---
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True, pady=5)
        
        # Tab 1: Standard Run
        tab_std = ttk.Frame(self.notebook)
        self.notebook.add(tab_std, text="Standard Run")
        self._setup_std_tab(tab_std)
        
        # Tab 2: PWM Ramp
        tab_ramp = ttk.Frame(self.notebook)
        self.notebook.add(tab_ramp, text="PWM Ramp")
        self._setup_ramp_tab(tab_ramp)
        
        # --- Log Output (Bottom) ---
        log_frame = ttk.LabelFrame(main_frame, text="Serial Log")
        log_frame.pack(fill="both", expand=True, pady=5)
        self.txt_log = scrolledtext.ScrolledText(log_frame, height=10, state=DISABLED)
        self.txt_log.pack(fill="both", expand=True)

    def _setup_std_tab(self, parent):
        f = ttk.Frame(parent, padding=10)
        f.pack(fill="both", expand=True)
        
        # Controls
        ttk.Button(f, text="Run Once (r)", command=lambda: self._send("r")).grid(row=0, column=0, pady=5)
        
        self.cont_var = StringVar(value="Stop")
        self.btn_cont = ttk.Button(f, text="Start Continuous Loop", command=self._toggle_continuous)
        self.btn_cont.grid(row=0, column=1, padx=10)
        
        # Status
        self.lbl_status = ttk.Label(f, text="Status: Idle")
        self.lbl_status.grid(row=1, column=0, columnspan=2, pady=10)

    def _setup_ramp_tab(self, parent):
        f = ttk.Frame(parent, padding=10)
        f.pack(fill="both", expand=True)
        
        # Settings
        cfg = ttk.Frame(f)
        cfg.pack(fill="x")
        
        ttk.Label(cfg, text="Channel:").pack(side="left")
        self.ramp_ch_var = StringVar(value="0")
        ttk.Entry(cfg, textvariable=self.ramp_ch_var, width=5).pack(side="left", padx=5)
        
        ttk.Label(cfg, text="Step Size:").pack(side="left")
        self.ramp_step_var = StringVar(value="5")
        ttk.Entry(cfg, textvariable=self.ramp_step_var, width=5).pack(side="left", padx=5)
        
        self.btn_ramp = ttk.Button(cfg, text="Start Ramp", command=self._start_ramp)
        self.btn_ramp.pack(side="left", padx=20)
        
        # Plot Area
        self.fig = Figure(figsize=(5, 3), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("PWM Ramp Response")
        self.ax.set_xlabel("PWM (0-255)")
        self.ax.set_ylabel("Voltage (V)")
        self.ax.grid(True)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=f)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True, pady=10)

    # --- Connectivity ---
    def _toggle_connect(self):
        if not self.is_connected:
            try:
                import serial
                self.ser = serial.Serial(self.port_var.get(), BAUD, timeout=READ_TIMEOUT)
                self.reader = SerialReader(self.ser, self.msg_queue)
                self.reader.start()
                self.is_connected = True
                self.btn_connect.config(text="Disconnect")
                self._log("Connected.")
            except Exception as e:
                messagebox.showerror("Error", str(e))
        else:
            if self.reader: self.reader.stop()
            if self.ser: self.ser.close()
            self.is_connected = False
            self.btn_connect.config(text="Connect")
            self._log("Disconnected.")

    def _send(self, cmd):
        if self.is_connected and self.ser:
            self.ser.write((cmd + "\n").encode())

    def _log(self, msg):
        self.txt_log.config(state=NORMAL)
        self.txt_log.insert(END, msg + "\n")
        self.txt_log.see(END)
        self.txt_log.config(state=DISABLED)

    # --- Ramp Logic ---
    def _start_ramp(self):
        if self.ramp_running: return
        try:
            ch = int(self.ramp_ch_var.get())
            step = int(self.ramp_step_var.get())
        except:
            messagebox.showerror("Input Error", "Invalid Channel or Step")
            return
            
        self.ramp_running = True
        self.ramp_data = []
        self.btn_ramp.config(state=DISABLED)
        self.ax.clear()
        self.ax.set_title(f"Ramp Channel {ch}")
        self.ax.grid(True)
        
        threading.Thread(target=self._ramp_worker, args=(ch, step), daemon=True).start()

    def _ramp_worker(self, ch, step):
        # 1. Dark Read
        self._send(f"pwm {ch} 0")
        time.sleep(0.1)
        
        pwm_vals = list(range(0, 256, step))
        if pwm_vals[-1] != 255: pwm_vals.append(255)
        
        for pwm in pwm_vals:
            if not self.is_connected: break
            
            # Set PWM
            self._send(f"pwm {ch} {pwm}")
            time.sleep(self.ramp_delay)
            
            # Request Readings
            self._send("a0")
            self._send("a1")
            
            # Wait for data (via queue in main loop)
            # This is tricky in threaded logic without blocking.
            # Simplified approach: We just send commands here, 
            # and the main loop parses and appends to data.
            # Ideally, we need a sync mechanism.
            
            # BETTER APPROACH: Do the logic here, read the queue here? 
            # No, queue is for GUI thread. 
            # Let's pause and assume the main loop catches the values.
            # We will use a shared "latest_readings" dict.
            
            time.sleep(0.05) # Wait for serial response
            
            # Grab latest from shared state (updated by gui_loop)
            a0 = self.latest_a0
            a1 = self.latest_a1
            
            self.ramp_data.append((pwm, a0, a1))
            
            # Trigger plot update
            self.root.after(0, self._update_plot)
            
        # Finish
        self._send(f"pwm {ch} 0") # Turn off
        self.ramp_running = False
        self.root.after(0, lambda: self.btn_ramp.config(state=NORMAL))
        self.root.after(0, self._save_ramp_csv)

    def _update_plot(self):
        if not self.ramp_data: return
        pwms = [x[0] for x in self.ramp_data]
        a1s = [x[2] for x in self.ramp_data]
        
        self.ax.clear()
        self.ax.plot(pwms, a1s, 'b-o')
        self.ax.set_xlabel("PWM")
        self.ax.set_ylabel("A1 Voltage")
        self.ax.grid(True)
        self.canvas.draw()

    def _save_ramp_csv(self):
        if not os.path.exists(SESSIONS_DIRNAME):
            os.makedirs(SESSIONS_DIRNAME)
        fname = f"{SESSIONS_DIRNAME}/Ramp_Ch{self.ramp_ch_var.get()}_{datetime.now().strftime('%H%M%S')}.csv"
        with open(fname, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["PWM", "A0", "A1"])
            writer.writerows(self.ramp_data)
        self._log(f"Ramp saved to {fname}")

    # --- Main Loop ---
    latest_a0 = 0.0
    latest_a1 = 0.0
    
    continuous_running = False
    
    def _toggle_continuous(self):
        self.continuous_running = not self.continuous_running
        if self.continuous_running:
            self.btn_cont.config(text="Stop Loop")
            threading.Thread(target=self._cont_worker, daemon=True).start()
        else:
            self.btn_cont.config(text="Start Continuous Loop")

    def _cont_worker(self):
        while self.continuous_running and self.is_connected:
            self._send("r")
            time.sleep(0.5) # Wait for run to finish (4 LEDs * 40ms * 2 ~ 320ms) + buffer

    def _gui_loop(self):
        while not self.msg_queue.empty():
            line = self.msg_queue.get()
            self._log(line)
            
            # Parse Values for Ramp
            # A0: 0.0450
            if "A0:" in line:
                m = PATTERNS["a0"].search(line)
                if m: self.latest_a0 = float(m.group(1))
            if "A1:" in line:
                m = PATTERNS["a1"].search(line)
                if m: self.latest_a1 = float(m.group(1))
                
        self.root.after(GUI_POLL_MS, self._gui_loop)

if __name__ == "__main__":
    root = Tk()
    app = App(root)
    root.mainloop()
