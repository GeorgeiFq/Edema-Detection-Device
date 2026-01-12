/***************************************************************
 * Photodiode Board PDV3 Control (With DAC offset cancel, DIFF_GAIN = 10)
 * - TLC59108 @ 0x40  (LED PWM)
 * - MAX30205 @ 0x48..0x4F (optional temperature)
 * - A0 = pre-diff (OPA2380 TIA output)
 * - A1 = post-diff (AD8629 diff-amp output)
 *
 * FIXES (2026-01-07):
 *  1) Per-channel DARK is now measured after forcing ALL LEDs OFF + settle.
 *     This prevents negative ΔA0/ΔA1 due to baseline contamination / recovery tails.
 *  2) vref=dark command now works (was previously unreachable).
 *  3) Serial input sanitization ignores non-printable bytes to avoid "Unknown command"
 *     caused by GUI junk bytes.
 ***************************************************************/
#include <Wire.h>
#include <math.h>

// =================== Version ===================
#define FW_VERSION "PDV3_FAST_TUNABLE_v10_FIXED"

// =================== User Config ===================
#define VREF_VOLTS 5.0f           // Arduino analog reference (UNO default 5.0 V)

// ADC averaging / timing (tunable at runtime)
static uint8_t  g_adc_avg_n         = 32;   // samples averaged per reported A0/A1
static uint16_t g_led_on_settle_ms  = 10;   // delay after LED turns on before sampling
static uint16_t g_led_off_settle_ms = 50;   // FIX: was 2 ms; default safer for recovery
static uint16_t g_dac_settle_ms     = 5;    // slightly increased for stability

// Auto Vref-from-dark settings
static float    g_auto_vref_factor  = 0.90f;  // Vref = factor * A0_dark (measured with DAC=0 V)
static uint16_t g_dark_settle_ms    = 100;    // FIX: default higher; you can still set via command

constexpr uint8_t NUM_CHANNELS = 4;

// I2C addresses
constexpr uint8_t TLC_ADDR = 0x40;   // TLC59108 LED driver
constexpr uint8_t DAC_ADDR = 0x4C;   // DAC7571 (A0=GND)
constexpr float   DAC_VDD  = 5.0f;   // DAC supply/reference

// ---- gslope config ----
static const uint8_t  GSLOPE_PWM_DEFAULTS[] = {120, 160, 200, 240};
static const uint8_t  GSLOPE_N = sizeof(GSLOPE_PWM_DEFAULTS) / sizeof(GSLOPE_PWM_DEFAULTS[0]);
static const uint16_t GSLOPE_SETTLE_MS = 150;   // settle after PWM change
static const float    A1_CEILING_V     = 4.5f;  // avoid railing during sweep

// ===== Auto-tune / chain constants =====
static const float DIFF_GAIN           = 10.0f;   // informational; used in safety check
static const float DELTA_A0_MIN_V      = 0.10f;
static const float DELTA_A0_MAX_V      = 0.20f;
static const uint8_t  PWM_MIN          = 5;
static const uint8_t  PWM_MAX          = 255;
static const uint8_t  TUNE_MAX_STEPS   = 7;
static const float    A1_MAX_SAFE_V    = 4.5f;

// =================== Forward Declarations ===================
static void     tlcSetPWM(uint8_t ch, uint8_t val);
static bool     tlcPresent();
static void     tlcInit();
static void     scanI2C();
static const char* nameForAddress(uint8_t addr);

static bool     dacPresent();
static void     dacWriteCode(uint16_t code);
static void     dacWriteVolts(float v);
static void     dacInit();

static bool     max30205_read_temp_at(uint8_t addr, float &outC);
static bool     max30205_probe();
static bool     maxPresent();
static float    readTemperatureC();

static float    readAnalogVoltage(uint8_t pin);
static inline float readA0();
static inline float readA1();

static void     autoVrefFromDarkA0();
static void     forceAllLedsOff();
static void     settleDarkBaseline();

static inline void ledOn(uint8_t ch, uint8_t pwm);
static inline void ledOff(uint8_t ch);

static float    tryPwmMeasureDeltaA0(uint8_t ch, uint8_t pwm,
                                     float &a0_dark, float &a1_dark, float &a1_on_out);
static uint8_t  autoTunePwmForDeltaA0(uint8_t ch, float &delta_a0,
                                      float &a0_dark_out, float &a1_dark_out, float &a1_on_out);

static void     autorun_basic();
static void     autorun_autotune();
static void     linfit(const float *x, const float *y, int n, float &m, float &b, float &r2);
static void     do_gslope_on_channel(uint8_t ch, const uint8_t *pwms, uint8_t np);

static void     printMenu();
static void     processCommand(String line);

static void     i2cPrintErr(uint8_t rc);

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
    case 0x40: return "TLC59108 LED driver (A3..A0=0000)";
    case 0x4C: return "DAC7571 (A0=GND)";
    case 0x4D: return "DAC7571 (A0=VDD)";
    case 0x48: return "MAX30205 Temp (A2=GND,A1=GND,A0=GND)";
    case 0x49: return "MAX30205 Temp";
    case 0x4A: return "MAX30205 Temp";
    case 0x4B: return "MAX30205 Temp (A2=GND,A1=VDD,A0=VDD)";
    case 0x4E: case 0x4F: return "MAX30205 Temp";
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
  i2cWrite8(TLC_ADDR, TLC_MODE1, 0x00);
  i2cWrite8(TLC_ADDR, TLC_MODE2, 0x00);
  i2cWrite8(TLC_ADDR, TLC_LEDOUT0, 0xAA); // OUT0..OUT3 -> individual PWM
  i2cWrite8(TLC_ADDR, TLC_LEDOUT1, 0x00);
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

// =================== DAC7571 ===================
static bool dacPresent() {
  Wire.beginTransmission(DAC_ADDR);
  return (Wire.endTransmission() == 0);
}

static void dacWriteCode(uint16_t code) {
  if (code > 0x0FFF) code = 0x0FFF;
  uint8_t msb = (uint8_t)((code >> 8) & 0x0F);
  uint8_t lsb = (uint8_t)(code & 0xFF);
  Wire.beginTransmission(DAC_ADDR);
  Wire.write(msb);
  Wire.write(lsb);
  Wire.endTransmission();
}

static void dacWriteVolts(float v) {
  if (v < 0) v = 0;
  if (v > DAC_VDD) v = DAC_VDD;
  uint16_t code = (uint16_t)lroundf((v / DAC_VDD) * 4095.0f);
  dacWriteCode(code);
}

static void dacInit() {
  dacWriteVolts(0.0f);
  delay(g_dac_settle_ms);
}

// =================== MAX30205 helpers ===================
static bool max30205_write_config(uint8_t addr, uint8_t cfg) {
  Wire.beginTransmission(addr);
  Wire.write((uint8_t)0x01);
  Wire.write(cfg);
  if (Wire.endTransmission(true) != 0) return false;
  return true;
}

static bool max30205_read_reg16(uint8_t addr, uint8_t reg, uint16_t &raw, bool useRepeatedStart) {
  Wire.beginTransmission(addr);
  Wire.write(reg);
  uint8_t rc = Wire.endTransmission(!useRepeatedStart);
  if (rc != 0) return false;

  uint8_t n = (uint8_t)Wire.requestFrom((uint8_t)addr, (uint8_t)2, (uint8_t)true);
  if (n != 2) return false;

  raw = ((uint16_t)Wire.read() << 8) | (uint16_t)Wire.read();
  return true;
}

static bool max30205_read_reg16_direct(uint8_t addr, uint16_t &raw) {
  uint8_t n = (uint8_t)Wire.requestFrom((uint8_t)addr, (uint8_t)2, (uint8_t)true);
  if (n != 2) return false;
  raw = ((uint16_t)Wire.read() << 8) | (uint16_t)Wire.read();
  return true;
}

static bool max30205_read_temp_at(uint8_t addr, float &outC) {
  uint16_t raw;
  if (max30205_read_reg16(addr, 0x00, raw, true)) {
    outC = ((int16_t)raw) / 256.0f;
    return true;
  }
  if (max30205_read_reg16(addr, 0x00, raw, false)) {
    outC = ((int16_t)raw) / 256.0f;
    return true;
  }
  if (max30205_read_reg16_direct(addr, raw)) {
    outC = ((int16_t)raw) / 256.0f;
    return true;
  }
  return false;
}

static bool max30205_probe() {
  const uint8_t preferred[] = {0x4B, 0x4A, 0x49, 0x48, 0x4C, 0x4D, 0x4E, 0x4F};
  for (uint8_t i = 0; i < sizeof(preferred); ++i) {
    uint8_t addr = preferred[i];
    Wire.beginTransmission(addr);
    if (Wire.endTransmission() == 0) {
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
  (void)analogRead(pin); // discard first conversion after mux switching
  uint32_t acc = 0;
  uint8_t n = g_adc_avg_n;
  if (n < 1) n = 1;
  if (n > 128) n = 128;
  for (uint8_t i = 0; i < n; i++) acc += (uint16_t)analogRead(pin);
  float avg = (float)acc / (float)n;
  return (avg / 1023.0f) * VREF_VOLTS;
}

static inline float readA0() { return readAnalogVoltage(A0); }
static inline float readA1() { return readAnalogVoltage(A1); }

// =================== LED helpers ===================
static inline void ledOn(uint8_t ch, uint8_t pwm) { tlcSetPWM(ch, pwm); }
static inline void ledOff(uint8_t ch)             { tlcSetPWM(ch, 0);   }

static void forceAllLedsOff() {
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ch++) ledOff(ch);
}

static void settleDarkBaseline() {
  // Allow analog chain to settle with all LEDs off
  delay(g_dark_settle_ms);
}

// =================== Auto Vref from dark (A0) ===================
static void autoVrefFromDarkA0() {
  forceAllLedsOff();

  // Force DAC to 0 V and settle
  dacWriteVolts(0.0f);
  delay(g_dac_settle_ms);

  settleDarkBaseline();

  float a0_dark = readA0();
  if (isnan(a0_dark) || a0_dark <= 0.0001f) {
    Serial.print(F("[WARN] AutoVref skipped: A0_dark invalid ("));
    if (isnan(a0_dark)) Serial.print(F("NaN"));
    else Serial.print(a0_dark, 4);
    Serial.println(F("). Leaving Vref=0 V."));
    return;
  }

  float vref = g_auto_vref_factor * a0_dark;
  if (vref < 0) vref = 0;
  if (vref > DAC_VDD) vref = DAC_VDD;

  dacWriteVolts(vref);
  delay(g_dac_settle_ms);

  // Give the diff stage a moment to re-center after the DAC step
  delay(10);

  float a1_dark_after = readA1();

  Serial.print(F("[AutoVref] A0_dark@0V=")); Serial.print(a0_dark, 4);
  Serial.print(F(" V, factor="));            Serial.print(g_auto_vref_factor, 3);
  Serial.print(F(" => Vref="));              Serial.print(vref, 4);
  Serial.print(F(" V, A1_dark(after)="));    Serial.print(a1_dark_after, 4);
  Serial.println(F(" V"));
}

// =================== Auto-tune helpers ===================
static float tryPwmMeasureDeltaA0(uint8_t ch, uint8_t pwm,
                                  float &a0_dark, float &a1_dark, float &a1_on_out) {
  // Ensure baseline is valid before each attempt
  forceAllLedsOff();
  delay(g_led_off_settle_ms);
  settleDarkBaseline();

  a0_dark = readA0();
  a1_dark = readA1();

  ledOn(ch, pwm);
  delay(g_led_on_settle_ms);

  float a0_on = readA0();
  float a1_on = readA1();

  ledOff(ch);
  delay(g_led_off_settle_ms);

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
      else                           return pwm;
    }

    if (hi < lo) break;
    pwm = (uint8_t)((lo + hi) / 2);
  }
  return best ? best : lo;
}

// =================== Linear fit (gslope) ===================
static void linfit(const float *x, const float *y, int n, float &m, float &b, float &r2) {
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

static void do_gslope_on_channel(uint8_t ch, const uint8_t *pwms, uint8_t np) {
  autoVrefFromDarkA0();

  Serial.println();
  Serial.print(F("=== gslope ch ")); Serial.print(ch); Serial.println(F(" ==="));

  // FIX: measure dark from a known, settled LED-off state
  forceAllLedsOff();
  delay(g_led_off_settle_ms);
  settleDarkBaseline();

  float a0_dark = readA0();
  float a1_dark = readA1();
  Serial.print(F("A0_dark: ")); Serial.println(a0_dark,4);
  Serial.print(F("A1_dark: ")); Serial.println(a1_dark,4);

  float x[10], y[10];
  uint8_t k = 0;

  for (uint8_t i=0; i<np; ++i) {
    uint8_t pwm = pwms[i];

    ledOn(ch, pwm);
    delay(GSLOPE_SETTLE_MS);
    float a0_on = readA0();
    float a1_on = readA1();
    ledOff(ch);
    delay(g_led_off_settle_ms);

    float dA0 = a0_on - a0_dark;
    float dA1 = a1_on - a1_dark;

    Serial.print(F("PWM ")); Serial.print(pwm);
    Serial.print(F(" -> A0_on: ")); Serial.print(a0_on,4);
    Serial.print(F(", A1_on: ")); Serial.print(a1_on,4);
    Serial.print(F(", dA0: "));   Serial.print(dA0,4);
    Serial.print(F(", dA1: "));   Serial.println(dA1,4);

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
  autoVrefFromDarkA0();

  Serial.println(F("=== Full Cycle Run (OUT0..OUT3 @ 255) ==="));

  float dv_a0[NUM_CHANNELS] = {0};
  float dv_a1[NUM_CHANNELS] = {0};

  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.print(F("\n=== Channel ")); Serial.print(ch); Serial.println(F(" ==="));

    // FIX: force known OFF state + settle BEFORE taking per-channel dark
    forceAllLedsOff();
    delay(g_led_off_settle_ms);
    settleDarkBaseline();

    float t_dark  = readTemperatureC();
    float v0_dark = readA0();
    float v1_dark = readA1();

    Serial.print(F("TEMP: "));
    if (isnan(t_dark)) Serial.println(F("N/A"));
    else { Serial.print(t_dark, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (dark): ")); Serial.println(v0_dark, 4);
    Serial.print(F("A1 (dark): ")); Serial.println(v1_dark, 4);

    // LED on
    ledOn(ch, 255);
    delay(g_led_on_settle_ms);

    float t_on  = readTemperatureC();
    float v0_on = readA0();
    float v1_on = readA1();

    Serial.print(F("TEMP (on): "));
    if (isnan(t_on)) Serial.println(F("N/A"));
    else { Serial.print(t_on, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (on): ")); Serial.println(v0_on, 4);
    Serial.print(F("A1 (on): ")); Serial.println(v1_on, 4);

    // LED off and settle
    ledOff(ch);
    delay(g_led_off_settle_ms);

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
  autoVrefFromDarkA0();

  Serial.println(F("=== Auto-tuned Autorun (ΔA0 window, with DAC offset cancel) ==="));

  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.println();
    Serial.print(F("=== Channel ")); Serial.print(ch); Serial.println(F(" (auto-tune) ==="));

    float dA0=0, a0_dark=0, a1_dark=0, a1_on_est=0;
    uint8_t tuned_pwm = autoTunePwmForDeltaA0(ch, dA0, a0_dark, a1_dark, a1_on_est);

    Serial.print(F("A0 (dark): ")); Serial.println(a0_dark, 4);
    Serial.print(F("A1 (dark): ")); Serial.println(a1_dark, 4);
    Serial.print(F("Tuned PWM: ")); Serial.println(tuned_pwm);
    Serial.print(F("ΔA0 (pre-check): ")); Serial.println(dA0, 4);

    ledOn(ch, tuned_pwm);
    delay(g_led_on_settle_ms);
    float a0_on = readA0();
    float a1_on = readA1();
    ledOff(ch);
    delay(g_led_off_settle_ms);

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
  Serial.println();
  Serial.println(F("Menu:"));
  Serial.println(F("1. Full Cycle Run          -> r"));
  Serial.println(F("2. Auto-tuned ΔA0 run      -> r0t"));
  Serial.println(F("3. Temperature reading     -> t"));
  Serial.println(F("4. Analog 0 reading        -> a0"));
  Serial.println(F("5. Analog 1 reading        -> a1"));
  Serial.println(F("Extras: i2c scan           -> scan"));
  Serial.println(F("        set PWM            -> pwm <ch 0..3> <0..255>"));
  Serial.println(F("        quick on/off       -> '<ch> on' or '<ch> off'"));
  Serial.println(F("        gslope             -> gslope [ch] [pwm1,pwm2,...]"));
  Serial.println(F("        set Vref (DAC)     -> vref <volts>   OR   vref=dark"));
  Serial.println(F("        AutoVref settings  -> vrefauto <factor>, darksettle <ms>"));
  Serial.println(F("        speed tuning       -> settle <ms>, avg <N>, dacsettle <ms>, offsettle <ms>"));
  Serial.println(F("        version            -> ver"));
  Serial.println(F("        help               -> h / help"));
  Serial.println();
}

static void i2cPrintErr(uint8_t rc) {
  switch (rc) {
    case 0: Serial.print(F("OK")); break;
    case 1: Serial.print(F("DATA_TOO_LONG")); break;
    case 2: Serial.print(F("NACK_ADDR")); break;
    case 3: Serial.print(F("NACK_DATA")); break;
    case 4: Serial.print(F("OTHER")); break;
    default: Serial.print(rc); break;
  }
}

static void processCommand(String line) {
  line.trim();
  line.toLowerCase();

  if (line.length() == 0) return;

  if (line == "ver") {
    Serial.print(F("FW: ")); Serial.println(FW_VERSION);
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

  // vref=dark (FIX: must be checked BEFORE startswith("vref"))
  if (line == "vref=dark" || line == "vref dark") {
    forceAllLedsOff();
    delay(g_led_off_settle_ms);
    settleDarkBaseline();

    float v = readA0();
    dacWriteVolts(v);
    delay(g_dac_settle_ms);

    Serial.print(F("Vref set to A0(dark) = ")); Serial.print(v,4); Serial.println(F(" V (DAC)"));
    return;
  }

  // vref <volts>
  if (line.startsWith("vref ")) {
    float v = line.substring(5).toFloat();
    dacWriteVolts(v);
    delay(g_dac_settle_ms);
    Serial.print(F("Vref set to ")); Serial.print(v,4); Serial.println(F(" V (DAC)"));
    return;
  }

  // vrefauto <factor>
  if (line.startsWith("vrefauto")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      float f = line.substring(s1+1).toFloat();
      if (f < 0.0f) f = 0.0f;
      if (f > 1.5f) f = 1.5f;
      g_auto_vref_factor = f;
      Serial.print(F("Auto Vref factor set to ")); Serial.println(g_auto_vref_factor, 3);
    } else {
      Serial.print(F("Auto Vref factor is ")); Serial.println(g_auto_vref_factor, 3);
      Serial.println(F("Usage: vrefauto <factor>"));
    }
    return;
  }

  // darksettle <ms>
  if (line.startsWith("darksettle")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      long ms = line.substring(s1+1).toInt();
      if (ms < 0) ms = 0;
      if (ms > 5000) ms = 5000;
      g_dark_settle_ms = (uint16_t)ms;
      Serial.print(F("Dark-settle set to ")); Serial.print(g_dark_settle_ms); Serial.println(F(" ms"));
    } else {
      Serial.print(F("Dark-settle is ")); Serial.print(g_dark_settle_ms); Serial.println(F(" ms"));
      Serial.println(F("Usage: darksettle <ms>"));
    }
    return;
  }

  // settle <ms> (LED-on settle)
  if (line.startsWith("settle")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      long ms = line.substring(s1+1).toInt();
      if (ms < 0) ms = 0;
      if (ms > 2000) ms = 2000;
      g_led_on_settle_ms = (uint16_t)ms;
      Serial.print(F("LED settle set to ")); Serial.print(g_led_on_settle_ms); Serial.println(F(" ms"));
    } else {
      Serial.print(F("LED settle is ")); Serial.print(g_led_on_settle_ms); Serial.println(F(" ms"));
      Serial.println(F("Usage: settle <ms>"));
    }
    return;
  }

  // offsettle <ms> (LED-off settle)
  if (line.startsWith("offsettle")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      long ms = line.substring(s1+1).toInt();
      if (ms < 0) ms = 0;
      if (ms > 5000) ms = 5000;
      g_led_off_settle_ms = (uint16_t)ms;
      Serial.print(F("LED off-settle set to ")); Serial.print(g_led_off_settle_ms); Serial.println(F(" ms"));
    } else {
      Serial.print(F("LED off-settle is ")); Serial.print(g_led_off_settle_ms); Serial.println(F(" ms"));
      Serial.println(F("Usage: offsettle <ms>"));
    }
    return;
  }

  // avg <N>
  if (line.startsWith("avg")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      long n = line.substring(s1+1).toInt();
      if (n < 1) n = 1;
      if (n > 128) n = 128;
      g_adc_avg_n = (uint8_t)n;
      Serial.print(F("ADC avg set to ")); Serial.print(g_adc_avg_n); Serial.println(F(" samples"));
    } else {
      Serial.print(F("ADC avg is ")); Serial.print(g_adc_avg_n); Serial.println(F(" samples"));
      Serial.println(F("Usage: avg <N>"));
    }
    return;
  }

  // dacsettle <ms>
  if (line.startsWith("dacsettle")) {
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      long ms = line.substring(s1+1).toInt();
      if (ms < 0) ms = 0;
      if (ms > 500) ms = 500;
      g_dac_settle_ms = (uint16_t)ms;
      Serial.print(F("DAC settle set to ")); Serial.print(g_dac_settle_ms); Serial.println(F(" ms"));
    } else {
      Serial.print(F("DAC settle is ")); Serial.print(g_dac_settle_ms); Serial.println(F(" ms"));
      Serial.println(F("Usage: dacsettle <ms>"));
    }
    return;
  }

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

  // gslope handler (kept from your original with minor cleanup)
  if (line.startsWith("gslope")) {
    uint8_t ch_start = 0, ch_end = NUM_CHANNELS-1;
    uint8_t pwms[10]; uint8_t np=0;

    int sp1 = line.indexOf(' ');
    if (sp1 > 0) {
      String arg1 = line.substring(sp1+1); arg1.trim();
      if (arg1.length()) {
        if (isDigit(arg1[0])) {
          int spaceOrComma = arg1.indexOf(' ');
          int commaOnly    = arg1.indexOf(',');
          if (spaceOrComma < 0 && commaOnly < 0) {
            int ch = arg1.toInt();
            if (ch >= 0 && ch < (int)NUM_CHANNELS) { ch_start=ch_end=(uint8_t)ch; }
          }
        }
        String list;
        int sp2 = arg1.indexOf(' ');
        if (sp2 >= 0) list = arg1.substring(sp2+1);
        else if (arg1.indexOf(',') >= 0) list = arg1;
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

  Serial.println(F("Unknown command. Type 'help' for menu."));
}

// =================== Setup / Loop ===================
void setup() {
  Serial.begin(115200);
  while (!Serial) { ; }
  Serial.println(F("\nPhotodiode Board PDV3 Control (With DAC offset cancel)\nFW: " FW_VERSION));

  Wire.begin();
  Wire.setClock(100000);

  if (max30205_probe()) {
    Serial.print(F("MAX30205 detected at 0x"));
    if (MAX_ADDR < 16) Serial.print('0');
    Serial.println(MAX_ADDR, HEX);
  } else {
    Serial.println(F("WARNING: MAX30205 not detected on 0x48–0x4F. Temp reading will be N/A."));
  }

  scanI2C();

  if (dacPresent()) {
    dacInit();
    Serial.println(F("DAC7571 detected and initialized (0x4C)."));
  } else {
    Serial.println(F("WARNING: DAC7571 not detected (0x4C). Offset cancel will not work."));
  }

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

    // FIX: ignore non-printable junk bytes that some GUIs/serial stacks can inject
    if (c == '\n' || c == '\r') {
      if (line.length()) { processCommand(line); line = ""; }
      continue;
    }
    if ((uint8_t)c < 32 || (uint8_t)c > 126) {
      // ignore control chars and extended bytes
      continue;
    }

    // handle backspace
    if (c == '\b') {
      if (line.length()) line.remove(line.length()-1);
      continue;
    }

    line += c;

    // guard line length
    if (line.length() > 80) {
      // flush if someone spams garbage without newlines
      line = "";
    }
  }
}
