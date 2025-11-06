/*
 * nRF5340 Edema Device — I2C Bus Scanner
 *
 * SDA = P1.02, SCL = P1.03  (bound to i2c2 via app.overlay)
 *
 * Build:
 *   west build -b nrf5340dk_nrf5340_cpuapp . --pristine
 * Flash:
 *   west flash
 */

#include <zephyr/kernel.h>
#include <zephyr/device.h>
#include <zephyr/devicetree.h>
#include <zephyr/drivers/i2c.h>
#include <zephyr/sys/printk.h>

#define SCAN_MIN_ADDR 0x03
#define SCAN_MAX_ADDR 0x77

/* Devicetree node labels (these must match your DTS/overlay) */
#define I2C0_NODE DT_NODELABEL(i2c0)
#define I2C1_NODE DT_NODELABEL(i2c1)
#define I2C2_NODE DT_NODELABEL(i2c2)

/* Forward declarations */
static const struct device *get_bus(const struct device *dev, const char *name);
static void scan_bus(const struct device *i2c, const char *tag);

/* ----------------------- helpers ----------------------- */
static const struct device *get_bus(const struct device *dev, const char *name)
{
    if (!device_is_ready(dev)) {
        printk("[ERROR] %s not ready\n", name);
        return NULL;
    }
    printk("[OK] %s ready\n", name);
    return dev;
}

static void scan_bus(const struct device *i2c, const char *tag)
{
    if (!i2c) {
        printk("[WARN] %s bus not available\n", tag);
        return;
    }

    printk("=== I2C Scan on %s ===\n", tag);
    uint8_t dummy;
    int found = 0;

    for (uint16_t addr = SCAN_MIN_ADDR; addr <= SCAN_MAX_ADDR; addr++) {
        int ret = i2c_read(i2c, &dummy, 1, addr); /* probe by read */
        if (ret == 0) {
            found++;
            const char *name = "";
            switch (addr) {
                case 0x40: name = " (TLC59108 LED Driver)"; break;
                case 0x48: name = " (MAX30205 / ADS1115)"; break;
                case 0x49: name = " (MAX30205 / ADS1115)"; break;
                case 0x4A: name = " (MAX30205 / ADS1115)"; break;
                case 0x4B: name = " (MAX30205 / ADS1115)"; break;
                case 0x4C: name = " (DAC7571)"; break;
                default: break;
            }
            printk("  Found device at 0x%02X%s\n", addr, name);
        }
    }

    if (found == 0) {
        printk("  [INFO] No I2C devices detected on %s.\n", tag);
    }
    printk("--- End of %s scan ---\n\n", tag);
}

/* ------------------------- main ------------------------ */
void main(void)
{
    printk("\n========================================\n");
    printk(" nRF5340 Edema Device — I2C Scanner Start\n");
    printk("========================================\n\n");

    const struct device *i2c0 = NULL;
    const struct device *i2c1 = NULL;
    const struct device *i2c2 = NULL;

#if DT_NODE_HAS_STATUS(I2C0_NODE, okay)
    i2c0 = get_bus(DEVICE_DT_GET(I2C0_NODE), "i2c0");
#endif
#if DT_NODE_HAS_STATUS(I2C1_NODE, okay)
    i2c1 = get_bus(DEVICE_DT_GET(I2C1_NODE), "i2c1");
#endif
#if DT_NODE_HAS_STATUS(I2C2_NODE, okay)
    i2c2 = get_bus(DEVICE_DT_GET(I2C2_NODE), "i2c2");
#endif

    k_sleep(K_MSEC(200)); /* let console come up */

    if (i2c0) scan_bus(i2c0, "i2c0");
    if (i2c1) scan_bus(i2c1, "i2c1");
    if (i2c2) scan_bus(i2c2, "i2c2");

    printk("Scan complete. Looping...\n");

    while (1) {
        k_sleep(K_SECONDS(5));
        if (i2c0) scan_bus(i2c0, "i2c0");
        if (i2c1) scan_bus(i2c1, "i2c1");
        if (i2c2) scan_bus(i2c2, "i2c2");
    }
}
