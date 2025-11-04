/***************************************************************
 * Photodiode Board v1 Control (No DAC, DIFF_GAIN = 12)
 * - TLC59108 @ 0x40  (LED PWM)
 * - MAX30205 @ 0x48..0x4F (optional temperature)
 * - A0 = pre-diff (OPA2380 TIA output)
 * - A1 = post-diff (AD8629 diff-amp output, ≈12× A0)
 *
 * Commands:
 *   r               -> autorun (OUT0..OUT3 once, PWM=255)
 *   r0t             -> autorun with ΔA0 auto-tune (per LED PWM)
 *   gslope [ch] [list] -> fit ΔA1 vs ΔA0 slope (gain) while avoiding rail
 *   t               -> temperature
 *   a0 / a1         -> analog reads
 *   scan            -> I2C scan
 *   pwm c v         -> set LED channel c=0..3 to value v=0..255
 *   "<c> on" / "<c> off" -> quick LED control
 ***************************************************************/
#include <Wire.h>
#include <math.h>   // fabs/fabsf/lround/isnan

//November 4th, 10x version

// =================== User Config ===================
#define VREF_VOLTS 5.0f           // Arduino analog reference (UNO default 5.0 V)
#define A_SAMPLES  16             // ADC oversamples per reading
#define STEP_DELAY_MS 250         // Delay used in basic autorun steps
constexpr uint8_t NUM_CHANNELS = 4;

// I2C addresses
constexpr uint8_t TLC_ADDR = 0x40;   // TLC59108 LED driver

// ---- gslope config ----
static const uint8_t  GSLOPE_PWM_DEFAULTS[] = {120, 160, 200, 240};
static const uint8_t  GSLOPE_N = sizeof(GSLOPE_PWM_DEFAULTS) / sizeof(GSLOPE_PWM_DEFAULTS[0]);
static const uint16_t GSLOPE_SETTLE_MS = 150;   // settle after PWM change
static const float    A1_CEILING_V     = 4.5f;  // avoid railing during sweep

// ===== Auto-tune / chain constants =====
static const float DIFF_GAIN           = 12.0f;   // AD8629 diff stage: 12k/1k => ~12×
static const float DELTA_A0_MIN_V      = 0.10f;   // target window for ΔA0
static const float DELTA_A0_MAX_V      = 0.20f;
static const uint16_t LED_ON_SETTLE_MS = 150;     // ≥5τ for 200k//0.1uF
static const uint8_t  PWM_MIN          = 5;
static const uint8_t  PWM_MAX          = 255;
static const uint8_t  TUNE_MAX_STEPS   = 7;       // "binary-ish" search depth
static const float    A1_MAX_SAFE_V    = 4.5f;    // keep post-diff away from rail

// =================== Forward Declarations (prototypes) ===================
// TLC / I2C
static void     tlcSetPWM(uint8_t ch, uint8_t val);
static bool     tlcPresent();
static void     tlcInit();
static void     scanI2C();
static const char* nameForAddress(uint8_t addr);

// MAX30205
static bool     max30205_read_temp_at(uint8_t addr, float &outC);
static bool     max30205_probe();
static bool     maxPresent();
static float    readTemperatureC();

// Analog
static float    readAnalogVoltage(uint8_t pin);
static inline float readA0();
static inline float readA1();

// LED convenience
static inline void ledOn(uint8_t ch, uint8_t pwm);
static inline void ledOff(uint8_t ch);

// Auto-tune
static float    tryPwmMeasureDeltaA0(uint8_t ch, uint8_t pwm,
                                     float &a0_dark, float &a1_dark, float &a1_on_out);
static uint8_t  autoTunePwmForDeltaA0(uint8_t ch, float &delta_a0,
                                      float &a0_dark_out, float &a1_dark_out, float &a1_on_out);

// Autoruns & gslope
static void     autorun_basic();
static void     autorun_autotune();
static void     linfit(const float *x, const float *y, int n, float &m, float &b, float &r2);
static void     do_gslope_on_channel(uint8_t ch, const uint8_t *pwms, uint8_t np);

// Shell
static void     printMenu();
static void     processCommand(String line);

// =================== TLC59108 regs ===================
constexpr uint8_t TLC_MODE1    = 0x00;
constexpr uint8_t TLC_MODE2    = 0x01;
constexpr uint8_t TLC_PWM0     = 0x02;  // PWM0..PWM7 at 0x02..0x09
constexpr uint8_t TLC_GRPPWM   = 0x0A;
constexpr uint8_t TLC_GRPFREQ  = 0x0B;
constexpr uint8_t TLC_LEDOUT0  = 0x0C;  // OUT0..OUT3 (2 bits each)
constexpr uint8_t TLC_LEDOUT1  = 0x0D;  // OUT4..OUT7 (2 bits each)

// =================== Globals ===================
static uint8_t MAX_ADDR = 0; // 0=unknown/not found

// =================== I2C helpers ===================
static inline void i2cWrite8(uint8_t dev, uint8_t reg, uint8_t val) {
  Wire.beginTransmission(dev);
  Wire.write(reg);
  Wire.write(val);
  Wire.endTransmission();
}

static const char* nameForAddress(uint8_t addr) {
  switch (addr) {
    case 0x40: return "TLC59108 LED driver";
    case 0x48: case 0x49: case 0x4A: case 0x4B:
    case 0x4C: case 0x4D: case 0x4E: case 0x4F:
      return "MAX30205 Temp (strap)";
    default:   return "Unknown";
  }
}

static void scanI2C() {
  Serial.println(F("\n--- I2C Scan ---"));
  bool any = false;
  for (uint8_t addr = 0x03; addr <= 0x77; ++addr) {
    Wire.beginTransmission(addr);
    if (Wire.endTransmission() == 0) {
      if (any) Serial.print(F(", "));
      any = true;
      Serial.print(F("0x"));
      if (addr < 16) Serial.print('0');
      Serial.print(addr, HEX);
      Serial.print(F(" ("));
      Serial.print(nameForAddress(addr));
      Serial.print(F(")"));
      delay(3);
    }
  }
  if (!any) Serial.print(F("No I2C devices found"));
  Serial.println();
  Serial.println(F("--- End Scan ---\n"));
}

// =================== TLC control ===================
static bool tlcPresent() {
  Wire.beginTransmission(TLC_ADDR);
  return (Wire.endTransmission() == 0);
}

static void tlcInit() {
  i2cWrite8(TLC_ADDR, TLC_MODE1, 0x00);   // MODE1 normal
  i2cWrite8(TLC_ADDR, TLC_MODE2, 0x00);   // MODE2 totem-pole
  i2cWrite8(TLC_ADDR, TLC_LEDOUT0, 0xAA); // OUT0..OUT3 -> individual PWM
  i2cWrite8(TLC_ADDR, TLC_LEDOUT1, 0x00); // OUT4..OUT7 off
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    i2cWrite8(TLC_ADDR, (uint8_t)(TLC_PWM0 + ch), 0x00);
  }
  i2cWrite8(TLC_ADDR, TLC_GRPPWM,  0xFF);
  i2cWrite8(TLC_ADDR, TLC_GRPFREQ, 0x00);
}

static void tlcSetPWM(uint8_t ch, uint8_t val) {
  if (ch >= NUM_CHANNELS) return;
  i2cWrite8(TLC_ADDR, (uint8_t)(TLC_PWM0 + ch), val);
}

// =================== MAX30205 (optional) ===================
static bool max30205_read_temp_at(uint8_t addr, float &outC) {
  Wire.beginTransmission(addr);
  Wire.write((uint8_t)0x00);              // temperature register
  if (Wire.endTransmission(false) != 0) return false;

  int n = Wire.requestFrom((int)addr, 2);
  if (n != 2) return false;

  uint16_t raw = ((uint16_t)Wire.read() << 8) | Wire.read();
  float tc = (int16_t)raw / 256.0f;       // 1 LSB = 1/256 C
  if (tc < -55.0f || tc > 150.0f) return false;
  outC = tc;
  return true;
}

static bool max30205_probe() {
  for (uint8_t addr = 0x48; addr <= 0x4F; ++addr) {
    float tc;
    if (max30205_read_temp_at(addr, tc)) {
      MAX_ADDR = addr;
      return true;
    }
  }
  MAX_ADDR = 0;
  return false;
}

static bool maxPresent() { return MAX_ADDR != 0; }

static float readTemperatureC() {
  if (!maxPresent()) return NAN;
  float tc;
  if (max30205_read_temp_at(MAX_ADDR, tc)) return tc;
  return NAN;
}

// =================== Analog helpers ===================
static float readAnalogVoltage(uint8_t pin) {
  float acc = 0;
  for (int i = 0; i < A_SAMPLES; i++) acc += analogRead(pin);
  acc /= (float)A_SAMPLES;
  return (acc / 1023.0f) * VREF_VOLTS;   // UNO 10-bit ADC
}

static inline float readA0() { return readAnalogVoltage(A0); } // pre-diff
static inline float readA1() { return readAnalogVoltage(A1); } // post-diff

// =================== LED helpers ===================
static inline void ledOn(uint8_t ch, uint8_t pwm) { tlcSetPWM(ch, pwm); }
static inline void ledOff(uint8_t ch)             { tlcSetPWM(ch, 0);   }

// =================== Auto-tune helpers ===================
static float tryPwmMeasureDeltaA0(uint8_t ch, uint8_t pwm,
                                  float &a0_dark, float &a1_dark, float &a1_on_out) {
  a0_dark = readA0();
  a1_dark = readA1();
  ledOn(ch, pwm);
  delay(LED_ON_SETTLE_MS);
  float a0_on = readA0();
  float a1_on = readA1();
  ledOff(ch);
  a1_on_out = a1_on;
  return (a0_on - a0_dark);
}

static uint8_t autoTunePwmForDeltaA0(uint8_t ch, float &delta_a0,
                                     float &a0_dark_out, float &a1_dark_out,
                                     float &a1_on_out) {
  uint8_t lo = PWM_MIN, hi = PWM_MAX, best = 0;
  float best_err = 1e9f;
  uint8_t pwm = (lo + hi) / 2;

  for (uint8_t step = 0; step < TUNE_MAX_STEPS; ++step) {
    float a0_dark, a1_dark, a1_on;
    float dA0 = tryPwmMeasureDeltaA0(ch, pwm, a0_dark, a1_dark, a1_on);

    if (a1_dark + DIFF_GAIN * dA0 > A1_MAX_SAFE_V) {
      hi = (pwm > PWM_MIN) ? (uint8_t)(pwm - 1) : PWM_MIN;
    } else {
      float mid = 0.5f * (DELTA_A0_MIN_V + DELTA_A0_MAX_V);
      float err = fabsf(dA0 - mid);
      if (err < best_err) {
        best_err = err; best = pwm;
        delta_a0 = dA0; a0_dark_out = a0_dark; a1_dark_out = a1_dark; a1_on_out = a1_on;
      }
      if (dA0 < DELTA_A0_MIN_V)      lo = (pwm < PWM_MAX) ? (uint8_t)(pwm + 1) : PWM_MAX;
      else if (dA0 > DELTA_A0_MAX_V) hi = (pwm > PWM_MIN) ? (uint8_t)(pwm - 1) : PWM_MIN;
      else                           return pwm; // in window
    }

    if (hi < lo) break;
    pwm = (uint8_t)((lo + hi) / 2);
  }
  return best ? best : lo;
}

// =================== Linear fit (gslope) ===================
static void linfit(const float *x, const float *y, int n, float &m, float &b, float &r2)
{
  float sx=0, sy=0, sxx=0, sxy=0, syy=0;
  for (int i=0;i<n;i++){ sx+=x[i]; sy+=y[i]; sxx+=x[i]*x[i]; sxy+=x[i]*y[i]; syy+=y[i]*y[i]; }
  float denom = n*sxx - sx*sx;
  if (fabs(denom) < 1e-9f) { m=0; b=0; r2=0; return; }
  m = (n*sxy - sx*sy) / denom;
  b = (sy - m*sx) / n;
  float ss_tot = syy - (sy*sy)/n;
  float ss_res = 0;
  for (int i=0;i<n;i++){ float yi = m*x[i] + b; float e = (y[i]-yi); ss_res += e*e; }
  r2 = (ss_tot <= 0) ? 1.0f : (1.0f - ss_res/ss_tot);
}

static void do_gslope_on_channel(uint8_t ch, const uint8_t *pwms, uint8_t np)
{
  Serial.println();
  Serial.print(F("=== gslope ch ")); Serial.print(ch); Serial.println(F(" ==="));

  float a0_dark = readA0();
  float a1_dark = readA1();
  Serial.print(F("A0_dark: ")); Serial.println(a0_dark,4);
  Serial.print(F("A1_dark: ")); Serial.println(a1_dark,4);

  float x[10], y[10]; // store ΔA0, ΔA1
  uint8_t k = 0;

  for (uint8_t i=0; i<np; ++i) {
    uint8_t pwm = pwms[i];

    // drive, read
    ledOn(ch, pwm);
    delay(GSLOPE_SETTLE_MS);
    float a0_on = readA0();
    float a1_on = readA1();
    ledOff(ch);

    float dA0 = a0_on - a0_dark;
    float dA1 = a1_on - a1_dark;

    Serial.print(F("PWM ")); Serial.print(pwm);
    Serial.print(F(" -> A0_on: ")); Serial.print(a0_on,4);
    Serial.print(F(", A1_on: ")); Serial.print(a1_on,4);
    Serial.print(F(", dA0: "));   Serial.print(dA0,4);
    Serial.print(F(", dA1: "));   Serial.println(dA1,4);

    // guardrails
    if (a1_on > A1_CEILING_V) { Serial.println(F("NOTE: A1 near ceiling; skipping")); continue; }
    if (fabs(dA0) < 0.003f)   { Serial.println(F("NOTE: dA0 too small; skipping")); continue; }

    x[k] = dA0; y[k] = dA1; k++;
  }

  if (k < 2) { Serial.println(F("Insufficient valid points for fit.")); return; }

  float m,b,r2;
  linfit(x,y,k,m,b,r2);
  Serial.print(F("gslope result — slope: ")); Serial.print(m,3);
  Serial.print(F(" , intercept: "));         Serial.print(b,3);
  Serial.print(F(" , R^2: "));               Serial.println(r2,4);
}

// =================== Autorun sequences ===================
static void autorun_basic() {
  Serial.println(F("=== Full Cycle Run (OUT0..OUT3 @ 255) ==="));

  float dv_a0[NUM_CHANNELS] = {0};
  float dv_a1[NUM_CHANNELS] = {0};

  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.print(F("\n=== Channel ")); Serial.print(ch); Serial.println(F(" ==="));

    float t_dark  = readTemperatureC();
    float v0_dark = readA0();
    float v1_dark = readA1();

    Serial.print(F("TEMP: "));
    if (isnan(t_dark)) Serial.println(F("N/A"));
    else { Serial.print(t_dark, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (dark): ")); Serial.println(v0_dark, 4);
    Serial.print(F("A1 (dark): ")); Serial.println(v1_dark, 4);

    // LED on @ full
    ledOn(ch, 255);
    delay(STEP_DELAY_MS);
    float t_on  = readTemperatureC();
    float v0_on = readA0();
    float v1_on = readA1();

    Serial.print(F("TEMP (on): "));
    if (isnan(t_on)) Serial.println(F("N/A"));
    else { Serial.print(t_on, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (on): ")); Serial.println(v0_on, 4);
    Serial.print(F("A1 (on): ")); Serial.println(v1_on, 4);

    // LED off, settle
    ledOff(ch);
    delay(STEP_DELAY_MS);
    float v0_post = readA0();
    float v1_post = readA1();

    float dV0 = v0_on - v0_dark;
    float dV1 = v1_on - v1_dark;

    if (fabs(dV0) > 1e-3) {
      float gain = dV1 / dV0;
      Serial.print(F("Gain (ΔA1/ΔA0): ")); Serial.println(gain, 3);
    } else {
      Serial.println(F("Gain (ΔA1/ΔA0): N/A (very small ΔA0)"));
    }

    dv_a0[ch] = dV0;
    dv_a1[ch] = dV1;

    Serial.print(F("A0 (post-off): ")); Serial.println(v0_post, 4);
    Serial.print(F("A1 (post-off): ")); Serial.println(v1_post, 4);

    Serial.print(F("ΔV = ")); Serial.print(dV0, 4); Serial.println(F(" V"));
  }

  // Summaries
  Serial.println(F("\n[SUMMARY] Full Cycle Run — AC Change per LED (A0: pre-diff amp)"));
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    long mv = lround(dv_a0[ch] * 1000.0f);
    Serial.print(F("LED ")); Serial.print(ch);
    Serial.print(F(" AC Change : ")); Serial.print(mv); Serial.println(F(" mV"));
  }
  Serial.println();

  Serial.println(F("[SUMMARY] Full Cycle Run — AC Change per LED (A1: post-diff amp)"));
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    long mv = lround(dv_a1[ch] * 1000.0f);
    Serial.print(F("LED ")); Serial.print(ch);
    Serial.print(F(" AC Change : ")); Serial.print(mv); Serial.println(F(" mV"));
  }
  Serial.println();

  Serial.println(F("[SUMMARY] Full Cycle Run — Gain per LED (ΔA1/ΔA0)"));
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.print(F("LED ")); Serial.print(ch); Serial.print(F(" Gain : "));
    if (fabs(dv_a0[ch]) > 1e-3) Serial.println(dv_a1[ch] / dv_a0[ch], 3);
    else                        Serial.println(F("N/A"));
  }
  Serial.println(F("\n=== Autorun complete ==="));
}

static void autorun_autotune() {
  Serial.println(F("=== Auto-tuned Autorun (ΔA0 window, no DAC) ==="));

  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.println();
    Serial.print(F("=== Channel ")); Serial.print(ch); Serial.println(F(" (auto-tune) ==="));

    float dA0=0, a0_dark=0, a1_dark=0, a1_on_est=0;
    uint8_t tuned_pwm = autoTunePwmForDeltaA0(ch, dA0, a0_dark, a1_dark, a1_on_est);

    Serial.print(F("A0 (dark): ")); Serial.println(a0_dark, 4);
    Serial.print(F("A1 (dark): ")); Serial.println(a1_dark, 4);
    Serial.print(F("Tuned PWM: ")); Serial.println(tuned_pwm);
    Serial.print(F("ΔA0 (pre-check): ")); Serial.println(dA0, 4);

    // Final measurement at tuned PWM
    ledOn(ch, tuned_pwm);
    delay(LED_ON_SETTLE_MS);
    float a0_on = readA0();
    float a1_on = readA1();
    ledOff(ch);

    float dV0 = a0_on - a0_dark;
    float dV1 = a1_on - a1_dark;

    Serial.print(F("A0 (on): "));  Serial.println(a0_on, 4);
    Serial.print(F("A1 (on): "));  Serial.println(a1_on, 4);
    Serial.print(F("ΔV0: "));      Serial.println(dV0, 4);
    Serial.print(F("ΔV1: "));      Serial.println(dV1, 4);

    if (fabs(dV0) > 1e-3) {
      float gain = dV1 / dV0;
      Serial.print(F("Gain (ΔA1/ΔA0): ")); Serial.println(gain, 3);
    } else {
      Serial.println(F("Gain (ΔA1/ΔA0): N/A (very small ΔA0)"));
    }
  }
  Serial.println(F("\n[INFO] Auto-tuned autorun complete."));
}

// =================== Shell ===================
static void printMenu() {
  Serial.println(F(
    "\nMenu:\n"
    "1. Full Cycle Run          -> r\n"
    "2. Auto-tuned ΔA0 run      -> r0t\n"
    "3. Temperature reading     -> t\n"
    "4. Analog 0 reading        -> a0\n"
    "5. Analog 1 reading        -> a1\n"
    "Extras: i2c scan           -> scan\n"
    "        help               -> h / help\n"
    "        set PWM            -> pwm <ch 0..3> <0..255>\n"
    "        quick on/off       -> '<ch> on' or '<ch> off'\n"
    "        gslope [ch] [pwm1,pwm2,...]\n"
  ));
}

static void processCommand(String line) {
  line.trim();
  line.toLowerCase();

  // --- gslope [ch] [pwm list] ---
  if (line.startsWith("gslope")) {
    uint8_t ch_start = 0, ch_end = NUM_CHANNELS-1;
    uint8_t pwms[10]; uint8_t np=0;

    int sp1 = line.indexOf(' ');
    if (sp1 > 0) {
      String arg1 = line.substring(sp1+1); arg1.trim();
      if (arg1.length()) {
        // try to parse channel
        if (isDigit(arg1[0])) {
          int spaceOrComma = arg1.indexOf(' ');
          int commaOnly    = arg1.indexOf(',');
          if (spaceOrComma < 0 && commaOnly < 0) {
            int ch = arg1.toInt();
            if (ch >= 0 && ch < (int)NUM_CHANNELS) { ch_start=ch_end=(uint8_t)ch; }
          }
        }
        // find CSV list (either second token or the token itself if it contains commas)
        String list;
        int sp2 = arg1.indexOf(' ');
        if (sp2 >= 0) {
          list = arg1.substring(sp2+1);
        } else if (arg1.indexOf(',') >= 0) {
          list = arg1;
        }
        list.trim();
        if (list.length()) {
          int pos=0;
          while (pos < list.length() && np < 10) {
            int comma = list.indexOf(',', pos);
            String tok = (comma<0) ? list.substring(pos) : list.substring(pos, comma);
            tok.trim();
            if (tok.length()) {
              int v = tok.toInt();
              if (v>=0 && v<=255) pwms[np++] = (uint8_t)v;
            }
            if (comma<0) break; else pos = comma+1;
          }
        }
      }
    }
    if (np==0) for (uint8_t i=0;i<GSLOPE_N;i++) pwms[np++] = GSLOPE_PWM_DEFAULTS[i];

    for (uint8_t ch=ch_start; ch<=ch_end; ++ch) do_gslope_on_channel(ch, pwms, np);
    return;
  }

  if (line == "r0t") { autorun_autotune(); return; }
  if (line == "r")   { autorun_basic();   return; }

  if (line == "t") {
    float tc = readTemperatureC();
    Serial.print(F("TEMP: "));
    if (isnan(tc)) Serial.println(F("N/A"));
    else { Serial.print(tc, 2); Serial.println(F(" C")); }
    return;
  }
  if (line == "a0") { float v=readA0(); Serial.print(F("A0: ")); Serial.print(v,4); Serial.println(F(" V")); return; }
  if (line == "a1") { float v=readA1(); Serial.print(F("A1: ")); Serial.print(v,4); Serial.println(F(" V")); return; }

  if (line == "scan") { scanI2C(); return; }
  if (line == "h" || line == "help" || line == "menu") { printMenu(); return; }

  // pwm <ch> <val>
  if (line.startsWith("pwm")) {
    int ch=-1, val=-1;
    int s1=line.indexOf(' ');
    if (s1>0) {
      int s2=line.indexOf(' ', s1+1);
      if (s2>s1) { ch=line.substring(s1+1,s2).toInt(); val=line.substring(s2+1).toInt(); }
    }
    if (ch>=0 && ch<(int)NUM_CHANNELS && val>=0 && val<=255) {
      tlcSetPWM((uint8_t)ch,(uint8_t)val);
      Serial.print(F("PWM set: ch=")); Serial.print(ch); Serial.print(F(" val=")); Serial.println(val);
    } else {
      Serial.println(F("Usage: pwm <ch 0..3> <0..255>"));
    }
    return;
  }

  // "<ch> on/off"
  int sp=line.indexOf(' ');
  if (sp>0) {
    String left=line.substring(0,sp), right=line.substring(sp+1);
    left.trim(); right.trim();
    if (right=="on" || right=="off") {
      int ch=left.toInt();
      if (ch>=0 && ch<(int)NUM_CHANNELS) {
        if (right=="on") { tlcSetPWM((uint8_t)ch,255); Serial.print(F("LED ")); Serial.print(ch); Serial.println(F(" -> ON (255)")); }
        else             { tlcSetPWM((uint8_t)ch,0);   Serial.print(F("LED ")); Serial.print(ch); Serial.println(F(" -> OFF (0)")); }
        return;
      }
    }
  }

  Serial.println(F("Unknown command. Type 'help' for menu."));
}

// =================== Setup / Loop ===================
void setup() {
  Serial.begin(115200);
  while (!Serial) { ; }
  Serial.println(F("\nPhotodiode Board v1 Control (No DAC, Gain=12x)"));

  Wire.begin();
  Wire.setClock(400000);

  if (max30205_probe()) {
    Serial.print(F("MAX30205 detected at 0x"));
    if (MAX_ADDR < 16) Serial.print('0');
    Serial.println(MAX_ADDR, HEX);
  } else {
    Serial.println(F("WARNING: MAX30205 not detected on 0x48–0x4F. Temp reading will be N/A."));
  }

  scanI2C();

  if (tlcPresent()) {
    tlcInit();
    Serial.println(F("TLC59108 initialized."));
  } else {
    Serial.println(F("WARNING: TLC59108 not detected (0x40). LED commands will do nothing."));
  }

  printMenu();
}

void loop() {
  static String line;
  while (Serial.available()) {
    char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      if (line.length()) { processCommand(line); line = ""; }
    } else {
      line += c;
    }
  }
}
