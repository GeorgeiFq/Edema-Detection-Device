#include <Wire.h>

// =================== User Config ===================
#define VREF_VOLTS 5.0f           // Change to 3.3f if your analog reference is 3.3 V
#define A_SAMPLES  16             // ADC averages for smoother readings
#define STEP_DELAY_MS 250         // 0.25 s per autorun step
constexpr uint8_t NUM_CHANNELS = 4; // We use OUT0..OUT3

// I2C known device addresses (update if your hardware changes)
constexpr uint8_t TLC_ADDR  = 0x40;  // TLC59108 LED driver
static uint8_t MAX_ADDR = 0;   // 0 = unknown until probed

// ===== Auto-ref settings =====
static bool   g_ref_auto = true;         // firmware will set DAC from A0_dark if true
static float  g_ref_margin_v = 0.010f;   // 10 mV below dark A0
static float  g_dac_floor_v  = 0.012f;   // avoid DAC zero-code neighborhood
static float  g_dac_max_v    = 4.90f;    // headroom below 5V
static float  g_dac_vdd      = 5.0f;     // DAC7571 supply (set 3.3f if on 3.3V)
static uint8_t g_dac_addr    = 0x4C;     // adjust if your DAC is at a different addr



// TLC59108 registers
constexpr uint8_t TLC_MODE1    = 0x00;
constexpr uint8_t TLC_MODE2    = 0x01;
constexpr uint8_t TLC_PWM0     = 0x02;  // PWM0..PWM7 live at 0x02..0x09
constexpr uint8_t TLC_GRPPWM   = 0x0A;
constexpr uint8_t TLC_GRPFREQ  = 0x0B;
constexpr uint8_t TLC_LEDOUT0  = 0x0C;  // OUT0..OUT3 (2 bits each)
constexpr uint8_t TLC_LEDOUT1  = 0x0D;  // OUT4..OUT7 (2 bits each)

// LEDOUT per-channel mode (2-bit):
// 0b00 = OFF, 0b01 = ON (no PWM), 0b10 = individual PWM, 0b11 = group PWM
enum LedOutMode : uint8_t { LED_OFF = 0, LED_ON = 1, LED_PWM = 2, LED_GRP = 3 };

// =================== I2C helpers ===================
static inline void i2cWrite8(uint8_t dev, uint8_t reg, uint8_t val) {
  Wire.beginTransmission(dev);
  Wire.write(reg);
  Wire.write(val);
  Wire.endTransmission();
}
constexpr uint8_t DAC_ADDR = 0x4C;  // adjust if your A0 pin changes

static bool dacPresent() {
  Wire.beginTransmission(DAC_ADDR);
  return (Wire.endTransmission() == 0);
}

static inline uint16_t _dac_code_from_volts(float v) {
  if (v < 0) v = 0;
  if (v > g_dac_vdd) v = g_dac_vdd;
  float codef = (v / g_dac_vdd) * 4095.0f;
  if (codef < 0) codef = 0;
  if (codef > 4095) codef = 4095;
  return (uint16_t)(codef + 0.5f);
}

static float readA0_avg(uint8_t n = A_SAMPLES) {
  uint32_t acc = 0;
  for (uint8_t i=0; i<n; ++i) {
    acc += analogRead(A0);
  }
  float adc = acc / (float)n;
  // If your ADC reference is VREF_VOLTS
  return (adc * VREF_VOLTS) / 1023.0f;   // Arduino UNO 10-bit
}


// Write DAC7571 (fast-mode I2C OK). Uses "fast mode, straight write" command.
static void dac_set_volts(float v) {
  if (v < g_dac_floor_v) v = g_dac_floor_v;
  if (v > g_dac_max_v)   v = g_dac_max_v;
  uint16_t code = _dac_code_from_volts(v);

  Wire.beginTransmission(g_dac_addr);
  // DAC7571 expects MSB then LSB of 12-bit code left-justified in two bytes:
  Wire.write((uint8_t)(code >> 4));            // D11..D4
  Wire.write((uint8_t)((code & 0xF) << 4));    // D3..D0 << 4
  Wire.endTransmission();

  // For your console/logs:
  Serial.print("DAC set to "); Serial.print(v, 3); Serial.println(" V");
}


// Try to read temperature at a specific address; returns true on success
static bool max30205_read_temp_at(uint8_t addr, float &outC) {
  // Point pointer to temperature register (0x00) and do a repeated start read
  Wire.beginTransmission(addr);
  Wire.write((uint8_t)0x00);
  if (Wire.endTransmission(false) != 0) return false;   // NACK on write

  int n = Wire.requestFrom((int)addr, 2);               // read 2 bytes
  if (n != 2) return false;

  uint16_t raw = ((uint16_t)Wire.read() << 8) | Wire.read();
  // MAX30205: 2's complement, 1 LSB = 1/256 C
  float tc = (int16_t)raw / 256.0f;
  // sanity range: -55..150 C (widen if you need)
  if (tc < -55.0f || tc > 150.0f) return false;

  outC = tc;
  return true;
}
// Writes a 12-bit DAC code (0–4095)
static void dacWriteRaw(uint16_t code) {
  if (code > 4095) code = 4095;
  uint8_t msb = (code >> 8) & 0x0F;
  uint8_t lsb = code & 0xFF;
  Wire.beginTransmission(DAC_ADDR);
  Wire.write(msb);  // upper 4 bits go here
  Wire.write(lsb);
  Wire.endTransmission();
}

// Convenience: set output voltage (0–VREF)
static void dacSetVoltage(float volts) {
  if (volts < 0) volts = 0;
  if (volts > VREF_VOLTS) volts = VREF_VOLTS;
  uint16_t code = (uint16_t)((volts / VREF_VOLTS) * 4095.0f + 0.5f);
  dacWriteRaw(code);
  Serial.print(F("DAC set to "));
  Serial.print(volts, 3);
  Serial.println(F(" V"));
}

// Try to read temperature from candidate addresses and pick the first that responds with a sane value
static int findMax30205() {
  for (uint8_t addr = 0x48; addr <= 0x4F; ++addr) {
    // point to temperature register 0x00
    Wire.beginTransmission(addr);
    Wire.write((uint8_t)0x00);
    if (Wire.endTransmission(false) != 0) continue;     // NACK

    // read two bytes
    if (Wire.requestFrom(addr, (uint8_t)2) != 2) continue;
    uint16_t raw = (Wire.read() << 8) | Wire.read();
    float tc = (int16_t)raw / 256.0f;

    // sanity check: human body-ish range; widen if needed
    if (tc > -55.0f && tc < 150.0f) {
      return addr;  // looks like a MAX30205 here
    }
  }
  return -1; // not found
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

static bool maxPresent() {
  return MAX_ADDR != 0;
}

static float readTemperatureC() {
  if (!maxPresent()) return NAN;
  float tc;
  if (max30205_read_temp_at(MAX_ADDR, tc)) return tc;
  return NAN;
}

static inline bool i2cRead1(uint8_t dev, uint8_t reg, uint8_t &out) {
  Wire.beginTransmission(dev);
  Wire.write(reg);
  if (Wire.endTransmission(false) != 0) return false;
  if (Wire.requestFrom(dev, (uint8_t)1) != 1) return false;
  out = Wire.read();
  return true;
}

// =================== TLC control ===================
static bool tlcPresent() {
  Wire.beginTransmission(TLC_ADDR);
  return (Wire.endTransmission() == 0);
}

static void tlcInit() {
  // MODE1: normal run
  i2cWrite8(TLC_ADDR, TLC_MODE1, 0x00);
  // MODE2: totem-pole
  i2cWrite8(TLC_ADDR, TLC_MODE2, 0x00);

  // Put OUT0..OUT3 into individual PWM mode (0b10 per channel → 0xAA)
  i2cWrite8(TLC_ADDR, TLC_LEDOUT0, 0xAA);
  // Leave OUT4..OUT7 unused/off
  i2cWrite8(TLC_ADDR, TLC_LEDOUT1, 0x00);

  // Initialize PWM to 0 for channels 0..3
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    i2cWrite8(TLC_ADDR, (uint8_t)(TLC_PWM0 + ch), 0x00);
  }
  // Optional group settings (unused here)
  i2cWrite8(TLC_ADDR, TLC_GRPPWM,  0xFF);
  i2cWrite8(TLC_ADDR, TLC_GRPFREQ, 0x00);
}

static void tlcSetPWM(uint8_t ch, uint8_t val) {
  if (ch >= NUM_CHANNELS) return;
  i2cWrite8(TLC_ADDR, (uint8_t)(TLC_PWM0 + ch), val);
}

// =================== MAX30205 ===================


// =================== Analog helpers ===================
static float readAnalogVoltage(uint8_t pin) {
  float acc = 0;
  for (int i = 0; i < A_SAMPLES; i++) acc += analogRead(pin);
  acc /= (float)A_SAMPLES;
  return (acc / 1023.0f) * VREF_VOLTS;
}

// =================== I2C scanner ===================
static const char* nameForAddress(uint8_t addr) {
  switch (addr) {
    case 0x40: return "TLC59108 LED driver";
    case 0x4C: return "DAC7571 DAC";
    case 0x48: case 0x49: case 0x4A: case 0x4B:
    case 0x4D: case 0x4E: case 0x4F:
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

// =================== Autorun sequence ===================
// =================== Autorun sequence ===================
// Replaces your existing autorun()
static void autorun() {
  Serial.println(F("=== Full Cycle Run (OUT0..OUT3 @ 255) ==="));

  // Store AC changes per channel for summaries (A0 = pre-diff, A1 = post-diff)
  float dv_a0[NUM_CHANNELS] = {0};
  float dv_a1[NUM_CHANNELS] = {0};

  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    Serial.print(F("\n=== Channel "));
    Serial.print(ch);
    Serial.println(F(" ==="));

    // --- DARK (LED OFF) ---
    float t_dark  = readTemperatureC();
    float v0_dark = readAnalogVoltage(A0);  // pre-diff
    float v1_dark = readAnalogVoltage(A1);  // post-diff

    Serial.print(F("TEMP: "));
    if (isnan(t_dark)) Serial.println(F("N/A"));
    else { Serial.print(t_dark, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (dark): ")); Serial.println(v0_dark, 4);
    Serial.print(F("A1 (dark): ")); Serial.println(v1_dark, 4);

    // --- ON (LED FULL) ---
    tlcSetPWM(ch, 255);
    delay(STEP_DELAY_MS);

    float t_on  = readTemperatureC();
    float v0_on = readAnalogVoltage(A0);    // pre-diff
    float v1_on = readAnalogVoltage(A1);    // post-diff

    Serial.print(F("TEMP (on): "));
    if (isnan(t_on)) Serial.println(F("N/A"));
    else { Serial.print(t_on, 2); Serial.println(F(" C")); }

    Serial.print(F("A0 (on): ")); Serial.println(v0_on, 4);
    Serial.print(F("A1 (on): ")); Serial.println(v1_on, 4);

    // --- POST (LED OFF) ---
    tlcSetPWM(ch, 0);
    delay(STEP_DELAY_MS);

    float v0_post = readAnalogVoltage(A0);
    float v1_post = readAnalogVoltage(A1);

    // Compute AC changes
    float dV0 = v0_on - v0_dark;  // A0 = pre-diff amp (your original ΔV)
    float dV1 = v1_on - v1_dark;  // A1 = post-diff amp

    // Save for summary
    dv_a0[ch] = dV0;
    dv_a1[ch] = dV1;

    Serial.print(F("A0 (post-off): ")); Serial.println(v0_post, 4);
    Serial.print(F("A1 (post-off): ")); Serial.println(v1_post, 4);

    // Keep your original ΔV line (based on A0)
    Serial.print(F("ΔV = "));
    Serial.print(dV0, 4);
    Serial.println(F(" V"));
  }

  // ---- Summaries ----
  // A0 (pre-diff amp)
  Serial.println(F("\n[SUMMARY] Full Cycle Run — AC Change per LED (A0: pre-diff amp)"));
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    long mv = lround(dv_a0[ch] * 1000.0f);
    Serial.print(F("LED "));
    Serial.print(ch);
    Serial.print(F(" AC Change : "));
    Serial.print(mv);
    Serial.println(F(" mV"));
  }
  Serial.println();

  // A1 (post-diff amp)
  Serial.println(F("[SUMMARY] Full Cycle Run — AC Change per LED (A1: post-diff amp)"));
  for (uint8_t ch = 0; ch < NUM_CHANNELS; ++ch) {
    long mv = lround(dv_a1[ch] * 1000.0f);
    Serial.print(F("LED "));
    Serial.print(ch);
    Serial.print(F(" AC Change : "));
    Serial.print(mv);
    Serial.println(F(" mV"));
  }
  Serial.println();

  Serial.println(F("=== Autorun complete ==="));
}



// =================== Command shell ===================
static void printMenu() {
  Serial.println(F(
    "\nMenu:\n"
    "1. Full Cycle Run          -> r\n"
    "2. LEDx Full brightness    -> 'x on' or 'x off'  (x=0..3)\n"
    "3. Temperature reading     -> t\n"
    "4. Analog 0 reading        -> a0\n"
    "5. Analog 1 reading        -> a1\n"
    "Extras: i2c scan           -> scan\n"
    "        help               -> h / help\n"
    "        set PWM            -> pwm <ch 0..3> <0..255>\n"   // <<< added
  ));
}


static void processCommand(String line) {
  line.trim();
  line.toLowerCase();

  if (line == "r") {
    autorun();
    return;
  }
  if (line == "t") {
    float tc = readTemperatureC();
    Serial.print(F("TEMP: "));
    if (isnan(tc)) Serial.println(F("N/A"));
    else { Serial.print(tc, 2); Serial.println(F(" C")); }
    return;
  }
  if (line == "a0") {
    float v = readAnalogVoltage(A0);
    Serial.print(F("A0: ")); Serial.print(v, 4); Serial.println(F(" V"));
    return;
  }
  if (line == "a1") {
    float v = readAnalogVoltage(A1);
    Serial.print(F("A1: ")); Serial.print(v, 4); Serial.println(F(" V"));
    return;
  }
  if (line == "scan") {
    scanI2C();
    return;
  }
  if (line == "h" || line == "help" || line == "menu") {
    printMenu();
    return;
  }
  // --- Per-channel PWM: "pwm <ch> <0..255>" ---
  if (line.startsWith("pwm")) {
    // expected forms: "pwm 0 128", "pwm 3 255", etc.
    int ch = -1, val = -1;

    // quick parse without tokenizers
    // find first space after "pwm"
    int s1 = line.indexOf(' ');
    if (s1 > 0) {
      int s2 = line.indexOf(' ', s1 + 1);
      if (s2 > s1) {
        ch  = line.substring(s1 + 1, s2).toInt();
        val = line.substring(s2 + 1).toInt();
      }
    }

    if (ch >= 0 && ch < (int)NUM_CHANNELS && val >= 0 && val <= 255) {
      tlcSetPWM((uint8_t)ch, (uint8_t)val);
      Serial.print(F("PWM set: ch=")); Serial.print(ch);
      Serial.print(F(" val=")); Serial.println(val);
    } else {
      Serial.println(F("Usage: pwm <ch 0..3> <0..255>"));
    }
    return;
  }

  // Parse "x on" / "x off"
  // Accept forms like: "0 on", "3 off"
  int sp = line.indexOf(' ');
  if (sp > 0) {
    String left = line.substring(0, sp);
    String right = line.substring(sp + 1);
    left.trim(); right.trim();
    if (right == "on" || right == "off") {
      int ch = left.toInt(); // if left isn’t numeric, toInt() returns 0; validate range below
      if (ch >= 0 && ch < (int)NUM_CHANNELS) {
        if (right == "on") {
          tlcSetPWM((uint8_t)ch, 255);
          Serial.print(F("LED ")); Serial.print(ch); Serial.println(F(" -> ON (255)"));
        } else {
          tlcSetPWM((uint8_t)ch, 0);
          Serial.print(F("LED ")); Serial.print(ch); Serial.println(F(" -> OFF (0)"));
        }
        return;
      }
    }
  }
  // --- DAC command ---
  if (line.startsWith("dac")) {
    // Example: "dac 2.5" or "dac 0x4C 2.5" or "dac raw 0x4C 2048"
    if (line.indexOf("raw") > 0) {
      // raw mode: dac raw 0x4C 2048
      int addrPos = line.indexOf("0x");
      int codePos = line.lastIndexOf(' ');
      if (addrPos > 0 && codePos > addrPos) {
        uint8_t addr = strtol(line.substring(addrPos, addrPos + 4).c_str(), NULL, 16);
        uint16_t code = line.substring(codePos + 1).toInt();
        Wire.beginTransmission(addr);
        Wire.write((code >> 8) & 0x0F);
        Wire.write(code & 0xFF);
        Wire.endTransmission();
        Serial.print(F("DAC raw "));
        Serial.print(line.substring(addrPos, addrPos + 4));
        Serial.print(F(" "));
        Serial.println(code);
      }
      return;
    }

    // Normal voltage mode
    float volts = 0;
    int space1 = line.indexOf(' ');
    if (space1 > 0) {
      String arg = line.substring(space1 + 1);
      arg.trim();

      // check if arg contains address
      if (arg.startsWith("0x")) {
        int nextSpace = arg.indexOf(' ', 0);
        if (nextSpace > 0) {
          String addrStr = arg.substring(0, nextSpace);
          uint8_t addr = strtol(addrStr.c_str(), NULL, 16);
          float v = arg.substring(nextSpace + 1).toFloat();
          uint16_t code = (uint16_t)((v / VREF_VOLTS) * 4095.0f + 0.5f);
          Wire.beginTransmission(addr);
          Wire.write((code >> 8) & 0x0F);
          Wire.write(code & 0xFF);
          Wire.endTransmission();
          Serial.print(F("DAC set to "));
          Serial.print(v, 3);
          Serial.println(F(" V"));
        }
      } else {
        volts = arg.toFloat();
        dacSetVoltage(volts);
      }
    } else {
      Serial.println(F("Usage: dac <volts>  or  dac 0x4C <volts>  or  dac raw 0x4C <code>"));
    }
    return;
  }
    // In processCommand(String line), before "Unknown command"
  if (line == "t?") {
    for (uint8_t addr = 0x48; addr <= 0x4F; ++addr) {
      float tc;
      bool ok = max30205_read_temp_at(addr, tc);
      Serial.print(F("Addr 0x"));
      if (addr < 16) Serial.print('0');
      Serial.print(addr, HEX);
      Serial.print(F(": "));
      if (ok) { Serial.print(tc, 2); Serial.println(F(" C")); }
      else    { Serial.println(F("N/A")); }
    }
    return;
  }


  Serial.println(F("Unknown command. Type 'help' for menu."));
}

// =================== Setup / Loop ===================
void setup() {
  Serial.begin(115200);
  while (!Serial) { ; }

  Serial.println(F("\nPhotodiode Board v1 Control"));
  Wire.begin();
  Wire.setClock(400000);
  // After Wire.begin(); Wire.setClock(400000);
  if (max30205_probe()) {
    Serial.print(F("MAX30205 detected at 0x"));
    if (MAX_ADDR < 16) Serial.print('0');
    Serial.println(MAX_ADDR, HEX);
  } else {
    Serial.println(F("WARNING: MAX30205 not detected on 0x48–0x4F. Temp reading will be N/A."));
  }


  // I2C scan at startup with friendly names
  scanI2C();

  // Initialize TLC if present
  if (tlcPresent()) {
    tlcInit();
    Serial.println(F("TLC59108 initialized."));
  } else {
    Serial.println(F("WARNING: TLC59108 not detected (0x40). LED commands will do nothing."));
  }

  if (maxPresent()) {
    Serial.println(F("MAX30205 detected."));
  } else {
    Serial.println(F("WARNING: MAX30205 not detected (0x4F). Temp reading will be N/A."));
  }

  if (dacPresent()) {
  Serial.println(F("DAC7571 detected."));
  } else {
    Serial.println(F("WARNING: DAC7571 not detected (0x4C)."));
  }
  if (dacPresent()) {
  Serial.println(F("DAC7571 detected."));
  } else {
    Serial.println(F("WARNING: DAC7571 not detected (0x4C)."));
  }



  printMenu();
}

void loop() {
  static String line;
  while (Serial.available()) {
    char c = (char)Serial.read();
    if (c == '\n' || c == '\r') {
      if (line.length()) {
        processCommand(line);
        line = "";
      }
    } else {
      line += c;
    }
  }
}
