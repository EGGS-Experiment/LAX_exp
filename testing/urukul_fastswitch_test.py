import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testidkurukul(EnvExperiment):
    """
    testidkurukul
    Testing urukul stuff
    """

    def build(self):
        self.frequency_mhz =                82.
        self.amplitude_pct =                50.

        self.phase_turns0 =                 0.
        # self.phase_turns1 =                 -0.325
        self.phase_turns1 =                 0.5

        self.time_delay_ns =                8

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

        # prepare device parameter values
        self._prepare_parameters()

    def _prepare_parameters(self):
        # calculate base waveform values
        self.ampl_asf0 =                    self.urukul0_ch0.amplitude_to_asf(self.amplitude_pct / 100.)
        self.ampl_asf1 =                    self.urukul0_ch0.amplitude_to_asf(self.amplitude_pct / 100. * 1.5)

        self.freq_ftw0 =                    self.urukul0_ch0.frequency_to_ftw(self.frequency_mhz * MHz)
        self.freq_ftw1 =                    self.urukul0_ch0.frequency_to_ftw(self.frequency_mhz * MHz)

        self.phase_ch0_offset_pow =         self.urukul0_ch0.turns_to_pow(self.phase_turns0)
        self.phase_ch1_offset_pow =         self.urukul0_ch0.turns_to_pow(self.phase_turns1)


        # calculate phase values to compensate for ch0 & ch1 delays
        self.phase_ch1_inherent_turns =     0.

        self.time_ch1_latency_ns =          1.66
        self.phase_ch1_latency_turns =      (self.frequency_mhz * MHz) * (self.time_ch1_latency_ns * ns)

        print(bin(np.int32(round(self.time_delay_ns))))
        time_delay_tmp_ns = np.int32(round(self.time_delay_ns)) & ~0x7
        print(bin(time_delay_tmp_ns))
        # self.phase_ch1_delay_turns =        (self.frequency_mhz * MHz) * (self.time_delay_ns * ns)
        self.phase_ch1_delay_turns =        (self.frequency_mhz * MHz) * (time_delay_tmp_ns * ns)
        # combine compensation values into ch1 phase: ch0 vs ch1 inherent phase shift/relation, ch0 vs ch1 inherent time delay, and ch1 offset
        self.phase_ch1_final_pow =          self.urukul0_ch0.turns_to_pow(self.phase_ch1_inherent_turns +
                                                                          self.phase_ch1_latency_turns +
                                                                          self.phase_ch1_delay_turns +
                                                                          self.phase_turns1)


        # calculate timings
        self.time_delay_mu =                self.core.seconds_to_mu(self.time_delay_ns * ns)
        self.time_phase_clear_mu =          self.core.seconds_to_mu(750 * ns)

        # tmp remove
        self.res00 = np.int32(0)
        self.res10 = np.int32(0)

        self.res01 = np.int32(0)
        self.res11 = np.int32(0)
        # tmp remove

    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        # prepare asynchronous/timing-agnostic stuff
        with parallel:
            # prepare - TTLs
            self.ttl8.off()
            self.ttl9.off()

            # prepare - DDSs
            self.urukul0_ch3.set_phase_mode(PHASE_MODE_ABSOLUTE)
            self.urukul1_ch3.set_phase_mode(PHASE_MODE_ABSOLUTE)

            self.urukul0_ch3.set_att(20 * dB)
            self.urukul1_ch3.set_att(10 * dB)

        # prepare - DDS profiles
        at_mu(now_mu() + 10000)
        with parallel:
            self.res00 = self.urukul0_ch3.set_mu(self.freq_ftw0, asf=0x01, pow_=0x0, profile=0)
            self.res10 = self.urukul1_ch3.set_mu(self.freq_ftw1, asf=0x01, pow_=self.phase_ch1_final_pow, profile=0)
        with parallel:
            self.res01 = self.urukul0_ch3.set_mu(self.freq_ftw0, asf=self.ampl_asf0, pow_=self.phase_ch0_offset_pow, profile=1)
            self.res11 = self.urukul1_ch3.set_mu(self.freq_ftw1, asf=self.ampl_asf1, pow_=self.phase_ch1_final_pow, profile=1)

        # set up registers
        at_mu(now_mu() + 10000)
        # self.set_cfr1_rdx(phase_autoclear=1, select_sine_output=1)
        with parallel:
            # self.urukul0_ch3.set_cfr1(phase_autoclear=1)
            # self.urukul1_ch3.set_cfr1(phase_autoclear=1)
            self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
            self.urukul1_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        with parallel:
            self.urukul0_ch3.set_cfr2(matched_latency_enable=1)
            self.urukul1_ch3.set_cfr2(matched_latency_enable=1)

        # set up for running
        at_mu(now_mu() + 10000)
        with parallel:
            self.urukul0_cpld.set_profile(0)
            self.urukul1_cpld.set_profile(0)
        with parallel:
            self.urukul0_cpld.cfg_switches(0b1000)
            self.urukul1_cpld.cfg_switches(0b1000)

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # prepare devices
        self._run_prepare()

        # prepare time
        # note: changing this delay causes relative phase shift; seems to follow according to frequency period
        # note: actually, seems to have more to do with SYNC_CLK/fixed timing periods
        delay_mu(125000)
        # note: we coarse align to previous SYNC_CLK period
        time_start_mu = now_mu() & ~7

        # RUN

        # start first waveform
        at_mu(time_start_mu)
        with parallel:
            self.urukul0_cpld.set_profile(1)
            # self.urukul1_cpld.set_profile(1)

            with sequential:
                delay_mu(475)
                self.ttl8.on()

                # tmp remove
                self.ttl9.on()
                # tmp remove

        at_mu(time_start_mu + self.time_delay_mu)
        delay_mu(475 + 4)
        self.ttl9.off()

        # start cancellation waveform
        # at_mu(time_start_mu + self.time_delay_mu)
        # # at_mu((time_start_mu + self.time_delay_mu) & ~0x3)
        # with parallel:
        #     self.urukul1_cpld.set_profile(1)
        #
        #     with sequential:
        #         delay_mu(475)
        #         self.ttl9.on()

        # at_mu(time_start_mu + 512)
        # with parallel:
        #     with sequential:
        #         self.urukul0_ch3.set_cfr1(phase_autoclear=1)
        #         self.urukul0_cpld.io_update.pulse_mu(8)
        #
        #     with sequential:
        #         delay_mu(787)
        #         self.ttl9.on()


        # # at_mu((time_start_mu + 508) & ~0x7)
        # at_mu((time_start_mu + 508))
        # with parallel:
        #     self.urukul0_cpld.io_update.pulse_mu(8)
        #
        #     with sequential:
        #         delay_mu(84)
        #         self.ttl9.on()



        # cleanup
        # with parallel:
        #     self.urukul0_cpld.set_profile(0)
        #     self.urukul1_cpld.set_profile(0)
        with parallel:
            self.urukul0_cpld.cfg_switches(0b0000)
            self.urukul1_cpld.cfg_switches(0b0000)
        delay_mu(100)
        with parallel:
            self.urukul0_cpld.set_profile(0)
            self.urukul1_cpld.set_profile(0)
            # self.urukul0_cpld.cfg_switches(0b0000)
            # self.urukul1_cpld.cfg_switches(0b0000)
            self.ttl8.off()
            self.ttl9.off()
        self.core.reset()


    @kernel(flags={"fast-math"})
    def set_cfr1_rdx(self,
                 power_down: TInt32 = 0b0000,
                 phase_autoclear: TInt32 = 0,
                 drg_load_lrr: TInt32 = 0,
                 drg_autoclear: TInt32 = 0,
                 phase_clear: TInt32 = 0,
                 internal_profile: TInt32 = 0,
                 ram_destination: TInt32 = 0,
                 ram_enable: TInt32 = 0,
                 manual_osk_external: TInt32 = 0,
                 osk_enable: TInt32 = 0,
                 select_auto_osk: TInt32 = 0,
                 select_sine_output: TInt32 = 1):
        """Set CFR1. See the AD9910 datasheet for parameter meanings.

        This method does not pulse IO_UPDATE.

        :param power_down: Power down bits.
        :param phase_autoclear: Autoclear phase accumulator.
        :param phase_clear: Asynchronous, static reset of the phase accumulator.
        :param drg_load_lrr: Load digital ramp generator LRR.
        :param drg_autoclear: Autoclear digital ramp generator.
        :param internal_profile: Internal profile control.
        :param ram_destination: RAM destination
            (:const:`RAM_DEST_FTW`, :const:`RAM_DEST_POW`,
            :const:`RAM_DEST_ASF`, :const:`RAM_DEST_POWASF`).
        :param ram_enable: RAM mode enable.
        :param manual_osk_external: Enable OSK pin control in manual OSK mode.
        :param osk_enable: Enable OSK mode.
        :param select_auto_osk: Select manual or automatic OSK mode.
        """
        with parallel:
            self.urukul0_ch3.write32(_AD9910_REG_CFR1,
                         (ram_enable << 31) |
                         (ram_destination << 29) |
                         (manual_osk_external << 23) |
                         (internal_profile << 17) |
                         (select_sine_output << 16) |
                         (drg_load_lrr << 15) |
                         (drg_autoclear << 14) |
                         (phase_autoclear << 13) |
                         (phase_clear << 11) |
                         (osk_enable << 9) |
                         (select_auto_osk << 8) |
                         (power_down << 4) |
                         2)  # SDIO input only, MSB first
            self.urukul1_ch3.write32(_AD9910_REG_CFR1,
                         (ram_enable << 31) |
                         (ram_destination << 29) |
                         (manual_osk_external << 23) |
                         (internal_profile << 17) |
                         (select_sine_output << 16) |
                         (drg_load_lrr << 15) |
                         (drg_autoclear << 14) |
                         (phase_autoclear << 13) |
                         (phase_clear << 11) |
                         (osk_enable << 9) |
                         (select_auto_osk << 8) |
                         (power_down << 4) |
                         2)  # SDIO input only, MSB first

    def analyze(self):
        print('\n\tresults:')
        print('\t\turukul0 profile 0 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res00)))
        print('\t\turukul1 profile 0 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res10)))
        print('\n')
        print('\t\turukul0 profile 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res01)))
        print('\t\turukul1 profile 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res11)))
        print('\n')

        print('\n\tio_update_delay:\t{}'.format(self.urukul0_ch3.sync_data.io_update_delay))
        print('\n\tsys_clk:\t{}'.format(self.urukul0_ch3.sysclk))
        print('\n\tsync_clk:\t{}'.format(self.urukul0_ch3.sysclk_per_mu * 4))
