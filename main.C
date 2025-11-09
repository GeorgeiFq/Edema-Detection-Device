/* SPDX-License-Identifier: LicenseRef-Nordic-5-Clause */
/* PD-2 menu + I2C scan + autorun over SEGGER RTT only (no UART/console) */

#include <zephyr/kernel.h>
#include <zephyr/device.h>
#include <zephyr/devicetree.h>
#include <zephyr/drivers/i2c.h>
#include <zephyr/drivers/adc.h>
#include <zephyr/logging/log.h>
#include <hal/nrf_saadc.h>

#include <SEGGER_RTT.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

LOG_MODULE_REGISTER(pd2_main, LOG_LEVEL_INF);

/* ===== I2C bus selection (alias first, then fallback) ===== */
#if DT_NODE_HAS_STATUS(DT_ALIAS(i2c1), okay)
#define I2C_NODE DT_ALIAS(i2c1)
#else
#define I2C_NODE DT_NODELABEL(arduino_i2c)
#endif

#ifndef CONFIG_I2C_SCAN_ADDR_START
#define CONFIG_I2C_SCAN_ADDR_START 0x03
#endif
#ifndef CONFIG_I2C_SCAN_ADDR_STOP
#define CONFIG_I2C_SCAN_ADDR_STOP  0x77
#endif

/* ===== ADC bind ===== */
#if DT_NODE_HAS_STATUS(DT_NODELABEL(adc), okay)
#define ADC_NODE DT_NODELABEL(adc)
#elif DT_NODE_HAS_STATUS(DT_NODELABEL(saadc), okay)
#define ADC_NODE DT_NODELABEL(saadc)
#else
#warning "No SAADC DT node labeled 'adc' or 'saadc'"
#define ADC_NODE DT_INVALID_NODE
#endif

/* Map your channels (AINx indices) */
#define ADC_CH_A0          0
#define ADC_CH_A1          1

#define ADC_RESOLUTION     12
#define ADC_GAIN           ADC_GAIN_1
#define ADC_REFERENCE      ADC_REF_INTERNAL
#define ACQ_TIME_SETTING   ADC_ACQ_TIME_DEFAULT

/* TLC59108 */
#define TLC_ADDR           0x40
#define TLC_REG_MODE1      0x00
#define TLC_REG_MODE2      0x01
#define TLC_REG_PWM0       0x02
#define TLC_REG_LEDOUT0    0x14
#define TLC_LEDOUT_PWM     0xAA

/* Autorun */
#define AUTORUN_PWM        255
#define AUTORUN_SETTLE_MS  250

/* Globals */
static const struct device *i2c_dev;
static const struct device *adc_dev;

/* ===== RTT tiny printf/IO helpers ======================================= */
static void rtt_puts(const char *s) { SEGGER_RTT_WriteString(0, s); }

static void rtt_printf(const char *fmt, ...)
{
    char buf[256];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    SEGGER_RTT_WriteString(0, buf);
}

/* Blocking line reader on RTT down channel 0 (CR/LF normalize, echoes) */
static void rtt_read_line(char *buf, size_t sz)
{
    size_t i = 0;
    for (;;) {
        /* poll for any bytes */
        char ch;
        int n = SEGGER_RTT_Read(0, &ch, 1);
        if (n == 0) { k_msleep(1); continue; }

        if (ch == '\r') ch = '\n';
        if (ch == '\n') { rtt_puts("\r\n"); break; }

        if (ch == 0x08 || ch == 0x7F) { /* backspace/delete */
            if (i > 0) { i--; rtt_puts("\b \b"); }
            continue;
        }

        if (isprint((unsigned char)ch)) {
            if (i + 1 < sz) { buf[i++] = ch; SEGGER_RTT_Write(0, &ch, 1); }
        }
    }
    buf[i] = '\0';
}

/* ===== I2C helpers ======================================================== */
static int probe_addr_quick(uint8_t addr)
{
    struct i2c_msg msg = { .buf = NULL, .len = 0, .flags = I2C_MSG_WRITE | I2C_MSG_STOP };
    return i2c_transfer(i2c_dev, &msg, 1, addr);
}

static int reg_read_u8(uint8_t addr, uint8_t reg, uint8_t *val)
{
    return i2c_write_read(i2c_dev, addr, &reg, 1, val, 1);
}
static int reg_read_n(uint8_t addr, uint8_t reg, uint8_t *buf, size_t n)
{
    return i2c_write_read(i2c_dev, addr, &reg, 1, buf, n);
}
static int reg_write_u8(uint8_t addr, uint8_t reg, uint8_t val)
{
    uint8_t tx[2] = {reg, val};
    return i2c_write(i2c_dev, tx, sizeof(tx), addr);
}

/* Classifiers */
static bool looks_like_TLC59108(uint8_t addr)
{
    uint8_t mode2=0, ledout0=0;
    if (reg_read_u8(addr, TLC_REG_MODE2, &mode2) != 0) return false;
    if (reg_read_u8(addr, TLC_REG_LEDOUT0, &ledout0) != 0) return false;
    return ((mode2 & 0x80u) == 0);
}

static bool looks_like_MAX30205(uint8_t addr)
{
    uint8_t cfg=0, t[2]={0};
    if (reg_read_u8(addr, 0x01, &cfg) != 0) return false;
    if (reg_read_n(addr, 0x00, t, 2) != 0)  return false;
    if ((cfg & 0xE0u) != 0) return false;
    if ((t[0]==0x00 && t[1]==0x00) || (t[0]==0xFF && t[1]==0xFF)) return false;
    return true;
}

static bool looks_like_ADS1x15(uint8_t addr)
{
    uint8_t cfg[2]={0}, conv[2]={0};
    if (reg_read_n(addr, 0x01, cfg, 2) != 0) return false;
    if (reg_read_n(addr, 0x00, conv, 2) != 0) return false;
    return !(cfg[0]==0xFF && cfg[1]==0xFF);
}

static const char *classify_addr(uint8_t addr)
{
    if (addr == TLC_ADDR && looks_like_TLC59108(addr)) return "TLC59108 LED driver";
    if (looks_like_MAX30205(addr))                      return "MAX30205 temp sensor";
    if (looks_like_ADS1x15(addr))                       return "ADS1x15 16 Bit ADC";
    return "Unknown (I2C responder)";
}

/* TLC59108 control */
static int tlc_init(void)
{
    if (reg_write_u8(TLC_ADDR, TLC_REG_MODE1, 0x00) != 0) return -EIO; /* wake */
    if (reg_write_u8(TLC_ADDR, TLC_REG_MODE2, 0x04) != 0) return -EIO; /* totem pole */
    if (reg_write_u8(TLC_ADDR, TLC_REG_LEDOUT0, TLC_LEDOUT_PWM) != 0) return -EIO;
    for (uint8_t ch=0; ch<4; ++ch) {
        if (reg_write_u8(TLC_ADDR, (uint8_t)(TLC_REG_PWM0 + ch), 0x00) != 0) return -EIO;
    }
    return 0;
}
static int tlc_set_pwm(uint8_t ch, uint8_t val)
{
    if (ch > 7) return -EINVAL;
    return reg_write_u8(TLC_ADDR, (uint8_t)(TLC_REG_PWM0 + ch), val);
}
static void tlc_all_off(void)
{
    for (uint8_t ch=0; ch<4; ++ch) (void)tlc_set_pwm(ch, 0);
}

/* MAX30205 read */
static int max30205_read_c(uint8_t addr, float *out_c)
{
    uint8_t t[2];
    if (reg_read_n(addr, 0x00, t, 2) != 0) return -EIO;
    int16_t raw = (int16_t)((t[0] << 8) | t[1]);
    *out_c = (float)raw / 256.0f;
    return 0;
}
static int read_any_max30205(float *out_c, uint8_t *addr_out)
{
    const uint8_t cands[] = {0x48,0x49,0x4A,0x4B,0x4C,0x4D,0x4E,0x4F};
    for (size_t i=0;i<sizeof(cands);++i) {
        if (probe_addr_quick(cands[i]) == 0 && looks_like_MAX30205(cands[i])) {
            if (max30205_read_c(cands[i], out_c) == 0) {
                if (addr_out) *addr_out = cands[i];
                return 0;
            }
        }
    }
    return -ENODEV;
}

/* ADC helpers */
static int adc_read_mv(uint8_t channel_id, int32_t *out_mv)
{
#if !DT_NODE_HAS_STATUS(ADC_NODE, okay)
    ARG_UNUSED(channel_id); ARG_UNUSED(out_mv);
    return -ENODEV;
#else
    if (!adc_dev) return -ENODEV;

    struct adc_channel_cfg cfg = {
        .gain             = ADC_GAIN,
        .reference        = ADC_REFERENCE,
        .acquisition_time = ACQ_TIME_SETTING,
#if defined(CONFIG_ADC_CONFIGURABLE_INPUTS)
        .channel_id       = channel_id,
        .input_positive   = NRF_SAADC_INPUT_AIN0 + channel_id,
#endif
        .differential     = 0,
    };

    int ret = adc_channel_setup(adc_dev, &cfg);
    if (ret) return ret;

    int16_t sample = 0;
    struct adc_sequence seq = {
        .channels    = BIT(channel_id),
        .buffer      = &sample,
        .buffer_size = sizeof(sample),
        .resolution  = ADC_RESOLUTION,
        .oversampling = 0,
        .calibrate = false,
    };

    ret = adc_read(adc_dev, &seq);
    if (ret) return ret;

    int32_t mv = sample;
    ret = adc_raw_to_millivolts(adc_ref_internal(adc_dev), ADC_GAIN, ADC_RESOLUTION, &mv);
    if (ret < 0) return ret;

    *out_mv = mv;
    return 0;
#endif
}

static float read_adc_or_nan(uint8_t ch, const char *label)
{
    int32_t mv = 0;
    int rc = adc_read_mv(ch, &mv);
    if (rc == 0) {
        return (float)mv / 1000.0f;
    } else {
        rtt_printf("%s N/A\r\n", label);
        return NAN;
    }
}

/* I2C scan */
static void do_scan(void)
{
    rtt_puts("\r\n===============================\r\n");
    rtt_puts("Scanning I2C bus...\r\n");
    rtt_puts("===============================\r\n");

    int found = 0;

    for (uint8_t addr = CONFIG_I2C_SCAN_ADDR_START; addr <= CONFIG_I2C_SCAN_ADDR_STOP; addr++) {
        if (probe_addr_quick(addr) == 0) {
            const char *name = classify_addr(addr);
            rtt_printf("Found 0x%02X  ->  %s\r\n", addr, name);
            found++;
        }
        k_msleep(2);
    }

    rtt_puts("\r\n===============================\r\n");
    rtt_printf("Scan complete. Devices found: %d\r\n", found);
    rtt_puts("Tips:\r\n");
    rtt_puts(" • TLC59108 usually at 0x40 on PD-2\r\n");
    rtt_puts(" • MAX30205 often at 0x48/0x49/0x4B (strapped)\r\n");
    rtt_puts(" • ADS1x15 ADCs vary (0x48–0x4B)\r\n");
    rtt_puts("===============================\r\n");
}

/* Autorun */
static void do_autorun(void)
{
    (void)tlc_init();

    rtt_puts("=== Full Cycle Run ===\r\n");

    for (uint8_t ch = 0; ch < 4; ++ch) {
        rtt_printf("=== Channel %u ===\r\n", ch);

        tlc_all_off();
        k_msleep(AUTORUN_SETTLE_MS);

        float a0_dark = read_adc_or_nan(ADC_CH_A0, "A0 (dark):");
        float a1_dark = read_adc_or_nan(ADC_CH_A1, "A1 (dark):");
        if (!isnan(a0_dark)) rtt_printf("A0 (dark): %.4f\r\n", a0_dark);
        if (!isnan(a1_dark)) rtt_printf("A1 (dark): %.4f\r\n", a1_dark);

        float temp_c = NAN;
        uint8_t taddr = 0;
        if (read_any_max30205(&temp_c, &taddr) == 0) rtt_printf("TEMP: %.2f\r\n", temp_c);
        else                                         rtt_puts("TEMP: N/A\r\n");

        (void)tlc_set_pwm(ch, AUTORUN_PWM);
        k_msleep(AUTORUN_SETTLE_MS);

        float a0_on = read_adc_or_nan(ADC_CH_A0, "A0 (on):");
        float a1_on = read_adc_or_nan(ADC_CH_A1, "A1 (on):");
        if (!isnan(a0_on)) rtt_printf("A0 (on): %.4f\r\n", a0_on);
        if (!isnan(a1_on)) rtt_printf("A1 (on): %.4f\r\n", a1_on);

        float temp_on_c = NAN;
        if (taddr && max30205_read_c(taddr, &temp_on_c) == 0)     rtt_printf("TEMP (on): %.2f\r\n", temp_on_c);
        else if (read_any_max30205(&temp_on_c, NULL) == 0)         rtt_printf("TEMP (on): %.2f\r\n", temp_on_c);
        else                                                       rtt_puts("TEMP (on): N/A\r\n");

        if (!isnan(a0_dark) && !isnan(a0_on)) {
            float dV = a0_on - a0_dark;
            rtt_printf("dV = %.6f\r\n", dV);
            if (!isnan(a1_dark) && !isnan(a1_on) && fabsf(dV) > 1e-6f) {
                float dA1 = a1_on - a1_dark;
                float gain = dA1 / dV;
                rtt_printf("Gain ( dA1/dA0 ) : %.6f\r\n", gain);
            }
        }

        (void)tlc_set_pwm(ch, 0);
        k_msleep(AUTORUN_SETTLE_MS);
        float a0_post = read_adc_or_nan(ADC_CH_A0, "A0 (post-off):");
        float a1_post = read_adc_or_nan(ADC_CH_A1, "A1 (post-off):");
        if (!isnan(a0_post)) rtt_printf("A0 (post-off): %.4f\r\n", a0_post);
        if (!isnan(a1_post)) rtt_printf("A1 (post-off): %.4f\r\n", a1_post);

        rtt_puts("\r\n");
    }

    rtt_puts("=== Autorun complete ===\r\n");
}

/* Menu */
static void print_menu_header(void)
{
    rtt_puts("\r\n--- I2C Devices Present (quick scan) ---\r\n");
    do_scan();
    rtt_puts("\r\nMenu:\r\n");
    rtt_puts("scan            -> scan\r\n");
    rtt_puts("Standard autorun-> r\r\n");
    rtt_puts("\r\nEnter command: ");
}

static void command_loop(void)
{
    char line[64];
    for (;;) {
        print_menu_header();
        rtt_read_line(line, sizeof(line));

        if (strcmp(line, "scan") == 0) {
            do_scan();
        } else if (strcmp(line, "r") == 0) {
            do_autorun();
        } else if (line[0] == '\0') {
            /* ignore */
        } else {
            rtt_printf("Unknown command: %s\r\n", line);
            rtt_puts("Valid: scan | r\r\n");
        }
    }
}

/* main */
void main(void)
{
    i2c_dev = DEVICE_DT_GET_OR_NULL(I2C_NODE);
    if (!i2c_dev || !device_is_ready(i2c_dev)) {
        LOG_ERR("I2C controller not ready.");
        return;
    }

#if DT_NODE_HAS_STATUS(ADC_NODE, okay)
    adc_dev = DEVICE_DT_GET_OR_NULL(ADC_NODE);
    if (adc_dev && !device_is_ready(adc_dev)) {
        LOG_WRN("ADC present but not ready; A0/A1 will be N/A.");
        adc_dev = NULL;
    }
#else
    adc_dev = NULL;
#endif

    rtt_printf("*** Booting nRF Connect SDK RTT-only build ***\r\n");
    rtt_printf("I2C controller: %s\r\n", i2c_dev->name);
#if DT_NODE_HAS_PROP(I2C_NODE, clock_frequency)
    rtt_printf("I2C bus clock: %u Hz\r\n", (uint32_t)DT_PROP(I2C_NODE, clock_frequency));
#endif

    k_msleep(200);
    command_loop();
}
