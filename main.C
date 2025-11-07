/* SPDX-License-Identifier: LicenseRef-Nordic-5-Clause */
/* I2C Scanner for Zephyr (RTT + printk hybrid) */

#include <zephyr/kernel.h>
#include <zephyr/device.h>
#include <zephyr/devicetree.h>
#include <zephyr/drivers/i2c.h>
#include <zephyr/logging/log.h>

LOG_MODULE_REGISTER(i2c_scanner, LOG_LEVEL_INF);

/* Pick I2C node from overlay alias or fallback */
#if DT_NODE_HAS_STATUS(DT_ALIAS(i2c1), okay)
#define I2C_NODE DT_ALIAS(i2c1)
#else
#define I2C_NODE DT_NODELABEL(arduino_i2c)
#endif

#ifndef CONFIG_I2C_SCAN_ADDR_START
#define CONFIG_I2C_SCAN_ADDR_START 0x03
#endif
#ifndef CONFIG_I2C_SCAN_ADDR_STOP
#define CONFIG_I2C_SCAN_ADDR_STOP 0x77
#endif

/* Zero-length write probe */
static int probe_addr_quick(const struct device *i2c, uint8_t addr)
{
    struct i2c_msg msg = {
        .buf = NULL,
        .len = 0,
        .flags = I2C_MSG_WRITE | I2C_MSG_STOP,
    };
    return i2c_transfer(i2c, &msg, 1, addr);
}

void main(void)
{
    /* Use OR_NULL so link never fails if the alias isn’t enabled */
    const struct device *i2c_dev = DEVICE_DT_GET_OR_NULL(I2C_NODE);

    LOG_INF("*** I2C Scanner Started ***");

    if (!i2c_dev || !device_is_ready(i2c_dev)) {
        LOG_ERR("I2C controller not ready (alias i2c1). Check app.overlay & CONFIG_I2C=y");
        return;
    }

    LOG_INF("Controller: %s", i2c_dev->name);

#if DT_NODE_HAS_PROP(I2C_NODE, clock_frequency)
    LOG_INF("Clock Frequency: %u Hz", (uint32_t)DT_PROP(I2C_NODE, clock_frequency));
#endif

    k_sleep(K_MSEC(200)); // allow bus settle

    printk("\n===============================\n");
    printk("Scanning I2C bus...\n");
    printk("===============================\n");

    int found = 0;

    for (uint8_t addr = CONFIG_I2C_SCAN_ADDR_START;
         addr <= CONFIG_I2C_SCAN_ADDR_STOP; addr++) {

        int ret = probe_addr_quick(i2c_dev, addr);
        if (ret == 0) {
            /* Print to RTT immediately (non-blocking relative to logging) */
            printk("Found device at 0x%02X\n", addr);
            found++;
        }

        /* Sleep briefly to avoid I2C saturation and allow RTT flush */
        k_sleep(K_MSEC(2));
    }

    printk("\n===============================\n");
    printk("Scan complete. Devices found: %d\n", found);
    printk("Reference: https://i2cdevices.org/addresses\n");
    printk("===============================\n");

    LOG_INF("I2C scan completed");
}
