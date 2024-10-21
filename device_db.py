
# Autogenerated for the ucla2 variant
core_addr = "192.168.1.76"
master_addr = "192.168.1.48"

device_db = {
    "core": {
        "type": "local",
        "module": "artiq.coredevice.core",
        "class": "Core",
        "arguments": {
            "host": core_addr,
            "ref_period": 1e-09,
            "analyzer_proxy": "core_analyzer",
            "target": "rv32g",
            "satellite_cpu_targets": {}
        },
    },
    "core_log": {
        "type": "controller",
        "host": "::1",
        # "host": master_addr,
        "port": 1068,
        "command": "aqctl_corelog -p {port} --bind {bind} " + core_addr
    },
    "core_moninj": {
        "type": "controller",
        # "host": "::1",
        "host": master_addr,
        "port_proxy": 1383,
        "port": 1384,
        "command": "aqctl_moninj_proxy --port-proxy {port_proxy} --port-control {port} --bind {bind} " + core_addr
    },
    "core_analyzer": {
        "type": "controller",
        # "host": "::1",
        "host": master_addr,
        "port_proxy": 1385,
        "port": 1386,
        "command": "aqctl_coreanalyzer_proxy --port-proxy {port_proxy} --port-control {port} --bind {bind} " + core_addr
    },
    "core_cache": {
        "type": "local",
        "module": "artiq.coredevice.cache",
        "class": "CoreCache"
    },
    "core_dma": {
        "type": "local",
        "module": "artiq.coredevice.dma",
        "class": "CoreDMA"
    },

    "i2c_switch0": {
        "type": "local",
        "module": "artiq.coredevice.i2c",
        "class": "I2CSwitch",
        "arguments": {"address": 0xe0}
    },
    "i2c_switch1": {
        "type": "local",
        "module": "artiq.coredevice.i2c",
        "class": "I2CSwitch",
        "arguments": {"address": 0xe2}
    },
}

# standalone peripherals

device_db["ttl0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000000},
}

device_db["ttl1"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000001},
}

device_db["ttl2"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000002},
}

device_db["ttl3"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000003},
}

device_db["ttl4"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000004},
}

device_db["ttl5"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000005},
}

device_db["ttl6"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000006},
}

device_db["ttl7"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLInOut",
    "arguments": {"channel": 0x000007},
}

device_db["ttl0_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x000008},
}

device_db["ttl1_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x000009},
}

device_db["ttl2_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000a},
}

device_db["ttl3_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000b},
}

device_db["ttl4_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000c},
}

device_db["ttl5_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000d},
}

device_db["ttl6_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000e},
}

device_db["ttl7_counter"] = {
    "type": "local",
    "module": "artiq.coredevice.edge_counter",
    "class": "EdgeCounter",
    "arguments": {"channel": 0x00000f},
}

device_db["ttl8"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000010},
}

device_db["ttl9"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000011},
}

device_db["ttl10"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000012},
}

device_db["ttl11"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000013},
}

device_db["ttl12"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000014},
}

device_db["ttl13"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000015},
}

device_db["ttl14"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000016},
}

device_db["ttl15"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000017},
}

device_db["eeprom_urukul0"] = {
    "type": "local",
    "module": "artiq.coredevice.kasli_i2c",
    "class": "KasliEEPROM",
    "arguments": {"port": "EEM2"}
}

device_db["spi_urukul0"] = {
    "type": "local",
    "module": "artiq.coredevice.spi2",
    "class": "SPIMaster",
    "arguments": {"channel": 0x000018}
}

device_db["ttl_urukul0_io_update"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000019}
}

device_db["ttl_urukul0_sw0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x00001a}
}

device_db["ttl_urukul0_sw1"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x00001b}
}

device_db["ttl_urukul0_sw2"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x00001c}
}

device_db["ttl_urukul0_sw3"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x00001d}
}

device_db["urukul0_cpld"] = {
    "type": "local",
    "module": "artiq.coredevice.urukul",
    "class": "CPLD",
    "arguments": {
        "spi_device": "spi_urukul0",
        "sync_device": None,
        "io_update_device": "ttl_urukul0_io_update",
        "refclk": 125000000.0,
        "clk_sel": 2,
        "clk_div": 0
    }
}

device_db["urukul0_ch0"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 4,
        "cpld_device": "urukul0_cpld",
        "sw_device": "ttl_urukul0_sw0"
    }
}

device_db["urukul0_ch1"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 5,
        "cpld_device": "urukul0_cpld",
        "sw_device": "ttl_urukul0_sw1"
    }
}

device_db["urukul0_ch2"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 6,
        "cpld_device": "urukul0_cpld",
        "sw_device": "ttl_urukul0_sw2"
    }
}

device_db["urukul0_ch3"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 7,
        "cpld_device": "urukul0_cpld",
        "sw_device": "ttl_urukul0_sw3"
    }
}

device_db["eeprom_urukul1"] = {
    "type": "local",
    "module": "artiq.coredevice.kasli_i2c",
    "class": "KasliEEPROM",
    "arguments": {"port": "EEM4"}
}

device_db["spi_urukul1"] = {
    "type": "local",
    "module": "artiq.coredevice.spi2",
    "class": "SPIMaster",
    "arguments": {"channel": 0x00001e}
}

device_db["ttl_urukul1_io_update"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x00001f}
}

device_db["ttl_urukul1_sw0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000020}
}

device_db["ttl_urukul1_sw1"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000021}
}

device_db["ttl_urukul1_sw2"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000022}
}

device_db["ttl_urukul1_sw3"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000023}
}

device_db["urukul1_cpld"] = {
    "type": "local",
    "module": "artiq.coredevice.urukul",
    "class": "CPLD",
    "arguments": {
        "spi_device": "spi_urukul1",
        "sync_device": None,
        "io_update_device": "ttl_urukul1_io_update",
        "refclk": 125000000.0,
        "clk_sel": 2,
        "clk_div": 0
    }
}

device_db["urukul1_ch0"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 4,
        "cpld_device": "urukul1_cpld",
        "sw_device": "ttl_urukul1_sw0"
    }
}

device_db["urukul1_ch1"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 5,
        "cpld_device": "urukul1_cpld",
        "sw_device": "ttl_urukul1_sw1"
    }
}

device_db["urukul1_ch2"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 6,
        "cpld_device": "urukul1_cpld",
        "sw_device": "ttl_urukul1_sw2"
    }
}

device_db["urukul1_ch3"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 7,
        "cpld_device": "urukul1_cpld",
        "sw_device": "ttl_urukul1_sw3"
    }
}

device_db["eeprom_urukul2"] = {
    "type": "local",
    "module": "artiq.coredevice.kasli_i2c",
    "class": "KasliEEPROM",
    "arguments": {"port": "EEM6"}
}

device_db["spi_urukul2"] = {
    "type": "local",
    "module": "artiq.coredevice.spi2",
    "class": "SPIMaster",
    "arguments": {"channel": 0x000024}
}

device_db["ttl_urukul2_io_update"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000025}
}

device_db["ttl_urukul2_sw0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000026}
}

device_db["ttl_urukul2_sw1"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000027}
}

device_db["ttl_urukul2_sw2"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000028}
}

device_db["ttl_urukul2_sw3"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000029}
}

device_db["urukul2_cpld"] = {
    "type": "local",
    "module": "artiq.coredevice.urukul",
    "class": "CPLD",
    "arguments": {
        "spi_device": "spi_urukul2",
        "sync_device": None,
        "io_update_device": "ttl_urukul2_io_update",
        "refclk": 125000000.0,
        "clk_sel": 2,
        "clk_div": 0
    }
}

device_db["urukul2_ch0"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 4,
        "cpld_device": "urukul2_cpld",
        "sw_device": "ttl_urukul2_sw0"
    }
}

device_db["urukul2_ch1"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 5,
        "cpld_device": "urukul2_cpld",
        "sw_device": "ttl_urukul2_sw1"
    }
}

device_db["urukul2_ch2"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 6,
        "cpld_device": "urukul2_cpld",
        "sw_device": "ttl_urukul2_sw2"
    }
}

device_db["urukul2_ch3"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_n": 32,
        "pll_en": 1,
        "chip_select": 7,
        "cpld_device": "urukul2_cpld",
        "sw_device": "ttl_urukul2_sw3"
    }
}

device_db["phaser0"] = {
    "type": "local",
    "module": "artiq.coredevice.phaser",
    "class": "Phaser",
    "arguments": {
        "channel_base": 0x00002a,
        "miso_delay": 1, "gw_rev": 1
    }
}

device_db["phaser1"] = {
    "type": "local",
    "module": "artiq.coredevice.phaser",
    "class": "Phaser",
    "arguments": {
        "channel_base": 0x00002f,
        "miso_delay": 1, "gw_rev": 1
    }
}

device_db["spi_sampler0_adc"] = {
    "type": "local",
    "module": "artiq.coredevice.spi2",
    "class": "SPIMaster",
    "arguments": {"channel": 0x000034}
}
device_db["spi_sampler0_pgia"] = {
    "type": "local",
    "module": "artiq.coredevice.spi2",
    "class": "SPIMaster",
    "arguments": {"channel": 0x000035}
}
device_db["ttl_sampler0_cnv"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000036},
}
device_db["sampler0"] = {
    "type": "local",
    "module": "artiq.coredevice.sampler",
    "class": "Sampler",
    "arguments": {
        "spi_adc_device": "spi_sampler0_adc",
        "spi_pgia_device": "spi_sampler0_pgia",
        "cnv_device": "ttl_sampler0_cnv",
        "hw_rev": "v2.2"
    }
}

device_db["led0"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000037}
}

device_db["led1"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000038}
}

device_db["led2"] = {
    "type": "local",
    "module": "artiq.coredevice.ttl",
    "class": "TTLOut",
    "arguments": {"channel": 0x000039}
}


##### STORE PHASER CONFIGURATION SETS #####

# original
# _phaser_trf_config = {
#     'lo_div_sel':           0b00,
#     'pll_div_sel':          0b01,
#     'rdiv':                 2,
#     'nint':                 23,
#     'prsc_sel':             0,
#     'cal_clk_sel':          0b1110
# }

# max
# _phaser_trf_config = {
#     'lo_div_sel':           0b11,
#     'pll_div_sel':          0b01,
#     'rdiv':                 21,
#     'nint':                 420,
#     'prsc_sel':             1,
#     'cal_clk_sel':          0b1010
# }


# customize phaser configuration
_phaser_dac_config = {
    'mixer_ena':            1,
    'nco_ena':              1,
    'fifo_offset':          5
}

# # TRF @ ~302.083853 MHz
_phaser_trf_config = {
    'pll_div_sel':          0b01,
    'rdiv':                 3,
    'nint':                 29,
    'prsc_sel':             0,
    'cal_clk_sel':          0b1101,

    'lo_div_sel':           0b11,
    'lo_div_bias':          0b00,
    'bufout_bias':          0b00,

    'tx_div_sel':           0b11,
    'tx_div_bias':          0b00
}


# # TRF @ 781.2512395 MHz
# _phaser_trf_config = {
#     'rdiv':                 2,
#     'nint':                 25,
#     'pll_div_sel':          0b01,
#     'prsc_sel':             0,
#
#     'icp':                  0b00000,
#     'icp_double':           0,
#
#     'cal_clk_sel':          0b1110,
#
#     'lo_div_sel':           0b11,
#     'lo_div_bias':          0b00,
#     'bufout_bias':          0b00,
#
#     'tx_div_sel':           0b10,
#     'tx_div_bias':          0b00
# }


# update phaser0 config options
device_db["phaser0"]["arguments"]["dac"] =      _phaser_dac_config
device_db["phaser0"]["arguments"]["trf0"] =     _phaser_trf_config
device_db["phaser0"]["arguments"]["trf1"] =     _phaser_trf_config

# update phaser1 config options
device_db["phaser1"]["arguments"]["dac"] =      _phaser_dac_config
device_db["phaser1"]["arguments"]["trf0"] =     _phaser_trf_config
device_db["phaser1"]["arguments"]["trf1"] =     _phaser_trf_config


