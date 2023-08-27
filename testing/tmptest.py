import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testidkurukul(EnvExperiment):
    """
    testidkurukul
    Testing
    """

    def build(self):
        self.frequency_mhz = 82.
        self.amplitude_pct = 50.

        self.phase_turns0 = 0.
        self.phase_turns1 = 0.5

        self.time_delay_ns = 20

    def prepare(self):
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # TTLs
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        # urukul CPLDs
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")

        # urukul channels
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")
        self.setattr_device("urukul2_ch2")
        self.setattr_device("urukul2_ch3")

        # prepare device values
        self._prepare_devices()

    def _prepare_devices(self):
        self.freq_ftw0 = self.urukul0_ch0.frequency_to_ftw(self.frequency_mhz * MHz)
        self.freq_ftw1 = self.urukul0_ch0.frequency_to_ftw(self.frequency_mhz * MHz)

        self.phase_pow0 = self.urukul0_ch0.turns_to_pow(self.phase_turns0)
        self.phase_pow1 = self.urukul0_ch0.turns_to_pow(self.phase_turns1)

        self.ampl_asf0 = self.urukul0_ch0.amplitude_to_asf(self.amplitude_pct / 100.)
        self.ampl_asf1 = self.urukul0_ch0.amplitude_to_asf(self.amplitude_pct / 100.)

        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_ns * ns)
        self.time_phase_clear_mu = self.core.seconds_to_mu(750 * ns)

    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        # prepare - TTLs
        self.ttl8.off()
        self.ttl9.off()

        # prepare - DDSs
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.urukul1_ch3.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.core.break_realtime()

        self.urukul0_ch3.set_att(14 * dB)
        self.urukul1_ch3.set_att(14 * dB)
        self.core.break_realtime()

        # prepare - profiles
        self.urukul0_ch3.set_mu(self.freq_ftw0, asf=0x01, pow_=0x0, profile=0)
        self.urukul1_ch3.set_mu(self.freq_ftw1, asf=0x01, pow_=0x0, profile=0)
        self.core.break_realtime()

        self.urukul0_ch3.set_mu(self.freq_ftw0, asf=self.ampl_asf0, pow_=self.phase_pow0, profile=1)
        self.urukul1_ch3.set_mu(self.freq_ftw1, asf=self.ampl_asf1, pow_=self.phase_pow1, profile=1)
        self.core.break_realtime()

        # set up registers
        self.urukul0_ch3.set_cfr1(phase_autoclear=1)
        self.urukul1_ch3.set_cfr1(phase_autoclear=1)

        self.urukul0_ch3.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch3.set_cfr2(matched_latency_enable=1)
        self.core.break_realtime()

        # set up for running
        self.urukul0_cpld.set_profile(0)
        self.urukul1_cpld.set_profile(0)
        self.core.break_realtime()

        self.urukul0_cpld.cfg_switches(0b1000)
        self.urukul1_cpld.cfg_switches(0b1000)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # prepare devices
        self._run_prepare()

        # prepare time
        delay_mu(1000000)
        time_start_mu = now_mu()

        # RUN
        # at_mu(time_start_mu)
        # with parallel:
        #     with sequential:
        #         self.urukul0_ch3.set_cfr1(phase_autoclear=1)
        #         self.urukul0_cpld.io_update.pulse_mu(8)
        #
        #     with sequential:
        #         self.urukul1_ch3.set_cfr1(phase_autoclear=1)
        #         self.urukul1_cpld.io_update.pulse_mu(8)

        at_mu(time_start_mu + self.time_phase_clear_mu)
        with parallel:
            self.urukul0_cpld.set_profile(1)

            with sequential:
                delay_mu(460)
                self.ttl8.on()


        at_mu(time_start_mu + self.time_phase_clear_mu + self.time_delay_mu)
        with parallel:
            self.ttl9.on()
            self.urukul1_cpld.set_profile(1)

        # cleanup
        # with parallel:
        #     self.urukul0_cpld.set_profile(0)
        #     self.urukul1_cpld.set_profile(0)
            # self.urukul0_cpld.cfg_switches(0b0000)
            # self.urukul1_cpld.cfg_switches(0b0000)
        delay_mu(100)
        with parallel:
            self.ttl8.off()
            self.ttl9.off()
        self.core.break_realtime()


    def analyze(self):
        pass
