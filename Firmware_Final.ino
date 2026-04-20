#include <Wire.h>

// ============================================================
// MPQ3326 Modular I2C Driver + Register Dump Utility
//
// CLI Commands:
//   help
//   dump
//   off
//   set <ch> <mA>        (ch=1..4, mA=0..50)  e.g. set 2 12.5
//   r<ch>  OR  r <ch>    (sets ch to 45mA, reads A0 50x avg, prints, turns ch off)
//   lin<ch> OR lin <ch>  (25..45mA step 1mA, each point reads A0 50x avg, prints table, turns ch off)
//   a0 [N]               (reads A0; default N=1 sample; prints counts + volts)
//   ao [N]               (alias for a0)
//   read [N]             (alias for a0)
//
// Notes:
// - ADC configured to 12-bit, assumed 3.3V reference.
// - Update ADC_REF_V if your effective AREF differs.
// ============================================================

// ---------- MPQ3326 register range ----------
static const uint8_t MPQ_REG_FIRST = 0x00;
static const uint8_t MPQ_REG_LAST  = 0x39;
static const uint8_t MPQ_REG_COUNT = (MPQ_REG_LAST - MPQ_REG_FIRST + 1);

// MPQ3326 possible 7-bit I2C addresses
static const uint8_t MPQ_ADDR_MIN = 0x30;
static const uint8_t MPQ_ADDR_MAX = 0x39;

// Set to 0 to auto-detect; otherwise hardcode your address (0x30..0x39)
static uint8_t g_mpqAddr = 0x30;

// Storage for register values (optional; used by dump)
static uint8_t g_regs[MPQ_REG_COUNT];

// ============================================================
// Current control constants (RISET-based)
// ============================================================
// You stated RISET = 10k to GND at ISET pin.
static const float RISET_KOHM = 10.0f;
// Datasheet: ILED_FS(mA) = 500 / RISET(kOhm)
static const float ISET_mA    = 500.0f / RISET_KOHM;  // -> 50 mA for 10k

// ============================================================
// A0 acquisition settings
// ============================================================
static const uint8_t  ANALOG_PIN        = A0;
static const uint16_t DEFAULT_A0_SAMPLES= 1;

// "rx" command settings
static const float    RX_CURRENT_mA     = 45.0f;
static const uint16_t RX_SAMPLES        = 10;

// "lin" command settings
static const int      LIN_START_mA      = 25;
static const int      LIN_END_mA        = 45;
static const int      LIN_STEP_mA       = 1;
static const uint16_t LIN_SAMPLES       = 10;

// Timing
static const uint16_t SETTLE_MS         = 10;  // wait after current set
static const uint16_t INTERSAMPLE_MS    = 2;   // delay between ADC samples

// ADC scaling (Nano 33 typically uses 3.3V ADC reference)
static const float    ADC_REF_V         = 3.3f;
static const uint8_t  ADC_BITS          = 12;
static const uint16_t ADC_MAX           = (1u << ADC_BITS) - 1u;

// ------------------------ Register labels ------------------------
const __FlashStringHelper* mpq3326RegLabel(uint8_t reg) {
  switch (reg) {
    case 0x00: return F("PWM Dimming Frequency Setting (FPWM[1:0])");
    case 0x01: return F("Control Register (FLTEN,LATCH,STH,SLEW,PS_EN,EN)");
    case 0x02: return F("Refresh + OTP Fault (FT_OTP, FRFSH[1:0])");
    case 0x03: return F("Refresh Frequency MSB (FRFSH[9:2])");
    case 0x04: return F("Channel Enable (CH9EN..CH16EN)");
    case 0x05: return F("Channel Enable (CH1EN..CH8EN)");
    case 0x06: return F("Open Fault Flags (CH9O..CH16O) [R]");
    case 0x07: return F("Open Fault Flags (CH1O..CH8O) [R]");
    case 0x08: return F("Short Fault Flags (CH9S..CH16S) [R]");
    case 0x09: return F("Short Fault Flags (CH1S..CH8S) [R]");

    case 0x0A: return F("ICH1 (Analog dim 6-bit)");
    case 0x0B: return F("PWM1 LSB (PWM1[3:0])");
    case 0x0C: return F("PWM1 MSB (PWM1[11:4])");

    case 0x0D: return F("ICH2 (Analog dim 6-bit)");
    case 0x0E: return F("PWM2 LSB (PWM2[3:0])");
    case 0x0F: return F("PWM2 MSB (PWM2[11:4])");

    case 0x10: return F("ICH3 (Analog dim 6-bit)");
    case 0x11: return F("PWM3 LSB (PWM3[3:0])");
    case 0x12: return F("PWM3 MSB (PWM3[11:4])");

    case 0x13: return F("ICH4 (Analog dim 6-bit)");
    case 0x14: return F("PWM4 LSB (PWM4[3:0])");
    case 0x15: return F("PWM4 MSB (PWM4[11:4])");

    default:   return F("Reserved/Unknown");
  }
}

// ============================================================
// Modular MPQ3326 I2C API
// ============================================================

bool mpq3326Ping(uint8_t addr) {
  Wire.beginTransmission(addr);
  return (Wire.endTransmission() == 0);
}

bool mpq3326WriteReg(uint8_t addr, uint8_t reg, uint8_t val) {
  Wire.beginTransmission(addr);
  Wire.write(reg);
  Wire.write(val);
  return (Wire.endTransmission() == 0);
}

bool mpq3326ReadReg(uint8_t addr, uint8_t reg, uint8_t &val) {
  Wire.beginTransmission(addr);
  Wire.write(reg);
  if (Wire.endTransmission(false) != 0) return false; // repeated start
  if (Wire.requestFrom((int)addr, 1) != 1) return false;
  val = Wire.read();
  return true;
}

bool mpq3326ReadAllRegs(uint8_t addr, uint8_t *outRegs, size_t outLen) {
  if (!outRegs || outLen < MPQ_REG_COUNT) return false;

  for (uint8_t reg = MPQ_REG_FIRST; reg <= MPQ_REG_LAST; reg++) {
    uint8_t v = 0;
    if (!mpq3326ReadReg(addr, reg, v)) return false;
    outRegs[reg - MPQ_REG_FIRST] = v;
  }
  return true;
}

bool mpq3326DumpAllRegs(uint8_t addr) {
  if (!mpq3326ReadAllRegs(addr, g_regs, MPQ_REG_COUNT)) {
    Serial.println(F("ERROR: mpq3326ReadAllRegs failed."));
    return false;
  }

  Serial.println(F("Full register dump (0x00..0x39):"));
  for (uint8_t reg = MPQ_REG_FIRST; reg <= MPQ_REG_LAST; reg++) {
    uint8_t val = g_regs[reg - MPQ_REG_FIRST];

    Serial.print(F("0x"));
    if (reg < 0x10) Serial.print('0');
    Serial.print(reg, HEX);

    Serial.print(F(" = 0x"));
    if (val < 0x10) Serial.print('0');
    Serial.print(val, HEX);

    Serial.print(F("  | "));
    Serial.println(mpq3326RegLabel(reg));
  }
  return true;
}

uint8_t mpq3326AutodetectAddr() {
  for (uint8_t a = MPQ_ADDR_MIN; a <= MPQ_ADDR_MAX; a++) {
    if (mpq3326Ping(a)) return a;
  }
  return 0;
}

// Writes 0x00 to registers 0x0A..0x15 (ICH1/PWM1..ICH4/PWM4)
bool mpq3326ZeroRegs_0A_to_15(uint8_t addr) {
  for (uint8_t reg = 0x0A; reg <= 0x15; reg++) {
    if (!mpq3326WriteReg(addr, reg, 0x00)) return false;
  }
  return true;
}

// ============================================================
// Channel helpers for analog/PWM dimming (CH1..CH4)
// ============================================================

static inline uint8_t mpqICHReg(uint8_t ch) {
  return (uint8_t)(0x0A + (ch - 1) * 3);
}
static inline uint8_t mpqPWMLsbReg(uint8_t ch) {
  return (uint8_t)(mpqICHReg(ch) + 1);
}
static inline uint8_t mpqPWMMsbReg(uint8_t ch) {
  return (uint8_t)(mpqICHReg(ch) + 2);
}

bool mpq3326SetPWM12(uint8_t addr, uint8_t ch, uint16_t pwm12) {
  if (ch < 1 || ch > 16) return false;
  if (pwm12 > 4095) pwm12 = 4095;

  uint8_t lsb = (uint8_t)(pwm12 & 0x0F);
  uint8_t msb = (uint8_t)((pwm12 >> 4) & 0xFF);

  // PWM updates when MSB is written => write LSB then MSB
  if (!mpq3326WriteReg(addr, mpqPWMLsbReg(ch), lsb)) return false;
  if (!mpq3326WriteReg(addr, mpqPWMMsbReg(ch), msb)) return false;
  return true;
}

bool mpq3326SetICH(uint8_t addr, uint8_t ch, uint8_t ichCode) {
  if (ch < 1 || ch > 16) return false;
  if (ichCode > 63) ichCode = 63;
  return mpq3326WriteReg(addr, mpqICHReg(ch), (uint8_t)(ichCode & 0x3F));
}

uint8_t mpq3326CurrentToICH(float mA) {
  if (mA <= 0.0f) return 0;
  if (mA >= ISET_mA) return 63;
  float codeF = 63.0f * (mA / ISET_mA);
  int code = (int)(codeF + 0.5f);
  if (code < 0) code = 0;
  if (code > 63) code = 63;
  return (uint8_t)code;
}

// Sets CH1..CH4 current in mA using analog dimming, PWM forced to 100% (4095)
bool mpq3326SetChannelCurrent_mA(uint8_t addr, uint8_t ch, float mA) {
  if (ch < 1 || ch > 4) return false;
  if (mA < 0.0f) mA = 0.0f;
  if (mA > ISET_mA) mA = ISET_mA;

  uint8_t ich = mpq3326CurrentToICH(mA);
  if (!mpq3326SetPWM12(addr, ch, 4095)) return false;
  if (!mpq3326SetICH(addr, ch, ich)) return false;
  return true;
}

bool mpq3326OffChannel(uint8_t addr, uint8_t ch) {
  if (ch < 1 || ch > 4) return false;
  if (!mpq3326SetPWM12(addr, ch, 0)) return false;
  if (!mpq3326SetICH(addr, ch, 0)) return false;
  return true;
}

bool mpq3326OffChannels1to4(uint8_t addr) {
  for (uint8_t ch = 1; ch <= 4; ch++) {
    if (!mpq3326OffChannel(addr, ch)) return false;
  }
  return true;
}

// ============================================================
// ADC helper: average a pin
// ============================================================
bool readAnalogAverage(uint8_t pin, uint16_t samples, float &avgCountsOut, float &avgVoltsOut) {
  if (samples < 1) samples = 1;
  uint32_t sum = 0;

  for (uint16_t i = 0; i < samples; i++) {
    int v = analogRead(pin);
    if (v < 0) v = 0;
    sum += (uint32_t)v;
    if (samples > 1) delay(INTERSAMPLE_MS);
  }

  avgCountsOut = (float)sum / (float)samples;
  avgVoltsOut  = (avgCountsOut * ADC_REF_V) / (float)ADC_MAX;
  return true;
}

// ============================================================
// "rx" helper
// ============================================================
bool runRx(uint8_t ch) {
  if (ch < 1 || ch > 4) return false;

  if (!mpq3326SetChannelCurrent_mA(g_mpqAddr, ch, RX_CURRENT_mA)) {
    Serial.println(F("ERR: rx set current failed"));
    mpq3326OffChannel(g_mpqAddr, ch);
    return false;
  }

  delay(SETTLE_MS);

  float avgCounts = 0.0f, avgVolts = 0.0f;
  readAnalogAverage(ANALOG_PIN, RX_SAMPLES, avgCounts, avgVolts);

  mpq3326OffChannel(g_mpqAddr, ch);

  Serial.print(F("RX CH"));
  Serial.print(ch);
  Serial.print(F(": I="));
  Serial.print(RX_CURRENT_mA, 2);
  Serial.print(F(" mA, A0_avg="));
  Serial.print(avgCounts, 2);
  Serial.print(F(" counts, "));
  Serial.print(avgVolts, 6);
  Serial.println(F(" V"));

  return true;
}

// ============================================================
// "lin" helper
// Output per line: "25mA CHx = y V"
// ============================================================
bool runLin(uint8_t ch) {
  if (ch < 1 || ch > 4) return false;

  for (int mA_int = LIN_START_mA; mA_int <= LIN_END_mA; mA_int += LIN_STEP_mA) {
    float mA = (float)mA_int;

    if (!mpq3326SetChannelCurrent_mA(g_mpqAddr, ch, mA)) {
      Serial.println(F("ERR: lin set current failed"));
      mpq3326OffChannel(g_mpqAddr, ch);
      return false;
    }

    delay(SETTLE_MS);

    float avgCounts = 0.0f, avgVolts = 0.0f;
    readAnalogAverage(ANALOG_PIN, LIN_SAMPLES, avgCounts, avgVolts);

    Serial.print(mA_int);
    Serial.print(F("mA CH"));
    Serial.print(ch);
    Serial.print(F(" = "));
    Serial.print(avgVolts, 6);
    Serial.println(F(" V"));
  }

  mpq3326OffChannel(g_mpqAddr, ch);
  return true;
}

// ============================================================
// Simple Serial CLI
// ============================================================

static char g_line[64];
static uint8_t g_lineLen = 0;

void printHelp() {
  Serial.println(F("Commands:"));
  Serial.println(F("  help"));
  Serial.println(F("  dump"));
  Serial.println(F("  off"));
  Serial.println(F("  set <ch> <mA>     (ch=1..4, mA=0..50) e.g. set 2 12.5"));
  Serial.println(F("  r<ch>  | r <ch>   (45mA, read A0 50x avg, print, off)"));
  Serial.println(F("  lin<ch>| lin <ch> (25..45mA step 1, read A0 50x avg, print table, off)"));
  Serial.println(F("  a0 [N]            (read A0; default N=1; prints counts + volts)"));
  Serial.println(F("  ao [N]            (alias for a0)"));
  Serial.println(F("  read [N]          (alias for a0)"));
  Serial.print(F("ISET full-scale (mA) = "));
  Serial.println(ISET_mA, 3);
}

static bool isAllDigits(const char *s) {
  if (!s || !*s) return false;
  while (*s) {
    if (*s < '0' || *s > '9') return false;
    s++;
  }
  return true;
}

void handleLine(const char *line) {
  char cmd[16] = {0};
  int a = 0;
  float b = 0.0f;

  // Generic parse: cmd + optional int + optional float
  int n = sscanf(line, "%15s %d %f", cmd, &a, &b);
  if (n < 1) {
    Serial.println(F("ERR: unknown command. Type 'help'."));
    return;
  }

  if (strcmp(cmd, "help") == 0) {
    printHelp();
    return;
  }

  if (strcmp(cmd, "dump") == 0) {
    mpq3326DumpAllRegs(g_mpqAddr);
    return;
  }

  if (strcmp(cmd, "off") == 0) {
    if (!mpq3326OffChannels1to4(g_mpqAddr)) {
      Serial.println(F("ERR: off failed"));
    } else {
      Serial.println(F("OK: channels 1-4 off"));
    }
    return;
  }

  // a0 / ao / read [N]
  if (strcmp(cmd, "a0") == 0 || strcmp(cmd, "ao") == 0 || strcmp(cmd, "read") == 0) {
    uint16_t samples = (n >= 2) ? (uint16_t)a : DEFAULT_A0_SAMPLES;
    if (samples < 1) samples = 1;

    float avgCounts = 0.0f, avgVolts = 0.0f;
    readAnalogAverage(ANALOG_PIN, samples, avgCounts, avgVolts);

    Serial.print(F("A0_avg("));
    Serial.print(samples);
    Serial.print(F(") = "));
    Serial.print(avgCounts, 2);
    Serial.print(F(" counts, "));
    Serial.print(avgVolts, 6);
    Serial.println(F(" V"));
    return;
  }

  // rX or r <ch>
  if (cmd[0] == 'r') {
    int rxCh = 0;

    if (strcmp(cmd, "r") == 0) {
      if (n < 2) {
        Serial.println(F("ERR: usage: r<ch> (e.g. r1) or r <ch>"));
        return;
      }
      rxCh = a;
    } else {
      const char *digits = &cmd[1];
      if (!isAllDigits(digits)) {
        Serial.println(F("ERR: usage: r<ch> (e.g. r1)"));
        return;
      }
      rxCh = atoi(digits);
    }

    if (rxCh < 1 || rxCh > 4) {
      Serial.println(F("ERR: channel must be 1..4"));
      return;
    }

    if (!runRx((uint8_t)rxCh)) {
      Serial.println(F("ERR: rx failed"));
    }
    return;
  }

  // linX or lin <ch>
  if (strncmp(cmd, "lin", 3) == 0) {
    int linCh = 0;

    if (strcmp(cmd, "lin") == 0) {
      if (n < 2) {
        Serial.println(F("ERR: usage: lin <ch> (e.g. lin 2) or lin<ch> (e.g. lin2)"));
        return;
      }
      linCh = a;
    } else {
      const char *digits = &cmd[3]; // after "lin"
      if (!isAllDigits(digits)) {
        Serial.println(F("ERR: usage: lin <ch> (e.g. lin 2) or lin<ch> (e.g. lin2)"));
        return;
      }
      linCh = atoi(digits);
    }

    if (linCh < 1 || linCh > 4) {
      Serial.println(F("ERR: channel must be 1..4"));
      return;
    }

    if (!runLin((uint8_t)linCh)) {
      Serial.println(F("ERR: lin failed"));
    }
    return;
  }

  // set <ch> <mA>
  if (strcmp(cmd, "set") == 0) {
    if (n < 3) {
      Serial.println(F("ERR: usage: set <ch> <mA>"));
      return;
    }
    int ch = a;
    float mA = b;

    if (ch < 1 || ch > 4) {
      Serial.println(F("ERR: channel must be 1..4"));
      return;
    }
    if (mA < 0.0f || mA > 50.0f) {
      Serial.println(F("ERR: mA must be 0..50"));
      return;
    }

    if (!mpq3326SetChannelCurrent_mA(g_mpqAddr, (uint8_t)ch, mA)) {
      Serial.println(F("ERR: set failed"));
      return;
    }

    uint8_t ich = mpq3326CurrentToICH(mA);
    Serial.print(F("OK: CH"));
    Serial.print(ch);
    Serial.print(F(" set to "));
    Serial.print(mA, 3);
    Serial.print(F(" mA (ICH="));
    Serial.print(ich);
    Serial.println(F(", PWM=4095)"));
    return;
  }

  Serial.println(F("ERR: unknown command. Type 'help'."));
}

void cliPoll() {
  while (Serial.available() > 0) {
    char c = (char)Serial.read();
    if (c == '\r') continue;

    if (c == '\n') {
      g_line[g_lineLen] = '\0';
      if (g_lineLen > 0) handleLine(g_line);
      g_lineLen = 0;
      return;
    }

    if (g_lineLen < (sizeof(g_line) - 1)) {
      g_line[g_lineLen++] = c;
    }
  }
}

// ============================================================
// Setup + Loop
// ============================================================

void setup() {
  Serial.begin(115200);
  while (!Serial) { delay(10); }

  analogReadResolution(ADC_BITS);

  Wire.begin();
  Wire.setClock(400000);

  Serial.println(F("MPQ3326 Modular Register Tool"));

  // Find device
  g_mpqAddr = (g_mpqAddr == 0) ? mpq3326AutodetectAddr() : g_mpqAddr;
  if (g_mpqAddr == 0) {
    Serial.println(F("ERROR: No MPQ3326 found at 0x30..0x39."));
    return;
  }

  Serial.print(F("Found MPQ3326 at 0x"));
  Serial.println(g_mpqAddr, HEX);

  Serial.print(F("Computed ISET full-scale (mA) from RISET=10k: "));
  Serial.println(ISET_mA, 3);

  mpq3326DumpAllRegs(g_mpqAddr);

  // Configure while EN=0
  if (!mpq3326WriteReg(g_mpqAddr, 0x01, 0x30)) {
    Serial.println(F("ERROR: Failed writing 0x30 to reg 0x01"));
    return;
  }

  // Disable CH9..CH16
  if (!mpq3326WriteReg(g_mpqAddr, 0x04, 0x00)) {
    Serial.println(F("ERROR: Failed writing 0x00 to reg 0x04"));
    return;
  }

  // Enable CH1..CH4, disable CH5..CH8
  if (!mpq3326WriteReg(g_mpqAddr, 0x05, 0x0F)) {
    Serial.println(F("ERROR: Failed writing 0x0F to reg 0x05"));
    return;
  }

  // Zero CH1..CH4 regs
  if (!mpq3326ZeroRegs_0A_to_15(g_mpqAddr)) {
    Serial.println(F("ERROR: Failed zeroing regs 0x0A..0x15"));
    return;
  }

  // Enable IC (EN=1)
  if (!mpq3326WriteReg(g_mpqAddr, 0x01, 0x31)) {
    Serial.println(F("ERROR: Failed writing 0x31 to reg 0x01"));
    return;
  }

  mpq3326DumpAllRegs(g_mpqAddr);

  Serial.println(F("Ready. Type 'help' for commands."));
  printHelp();
}

void loop() {
  cliPoll();
}
