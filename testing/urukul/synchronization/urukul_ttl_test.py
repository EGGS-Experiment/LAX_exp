import numpy as np
from artiq.experiment import *

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1


class UrukulTTLTest(EnvExperiment):
    """
    Urukul TTL Test
    Testing use of TTL to fast switch Urukul pulses.
    """

    def build(self):
        self.frequency_mhz =                350.
        self.amplitude_pct =                88.
        self.att_db =                       0.

        self.phase_turns0 =                 0.
        self.phase_turns1 =                 0.

        self.time_delay_ns =                1000

    def prepare(self):
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # TTLs
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        # urukul devices
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        # prepare device parameter values
        self._prepare_parameters()

    def _prepare_parameters(self):
        # calculate base waveform values
        self.ampl_asf0 =                    self.urukul0_ch3.amplitude_to_asf(self.amplitude_pct / 100.)
        self.freq_ftw0 =                    self.urukul0_ch3.frequency_to_ftw(self.frequency_mhz * MHz)

    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        # prepare asynchronous/timing-agnostic stuff
        # with parallel:
        # prepare - TTLs
        self.ttl8.off()
        self.ttl9.off()

        # prepare - DDSs
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_att(self.att_db * dB)

        # prepare - DDS profiles
        at_mu(now_mu() + 10000)
        self.urukul0_ch3.set_mu(self.freq_ftw0, asf=self.ampl_asf0, profile=0,
                                pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)

        # set up registers
        at_mu(now_mu() + 10000)
        self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch3.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch3.sw.off()


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()
        self._run_prepare()

        delay_mu(125000)
        at_mu(now_mu() & ~7)
        with parallel:
            self.urukul0_cpld.set_profile(0)
            self.urukul0_ch3.sw.on()

        delay_mu(1000)
        with parallel:
            self.ttl8.on()
            self.ttl9.on()

        delay_mu(8)
        with parallel:
            self.ttl8.off()
            self.ttl9.off()

        delay_mu(500)
        self.urukul0_ch3.sw.off()
        # self.core.reset()
