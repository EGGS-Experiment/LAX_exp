import numpy as np
from artiq.experiment import *

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1


class UrukulTrackingTest(EnvExperiment):
    """
    Urukul Tracking Test
    Testing phase tracking.
    """

    def build(self):
        self.frequency_mhz =                82.
        self.amplitude_pct =                50.
        self.att_db =                       10.

        self.phase_turns0 =                 0.
        self.phase_turns1 =                 0.5

        self.time_delay_ns =                8


    def prepare(self):
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # TTLs
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        # urukul devices
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

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
            self.urukul1_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
            self.urukul1_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)

            self.urukul1_ch1.set_att(self.att_db * dB)
            self.urukul1_ch2.set_att(self.att_db * dB)

        # prepare - DDS profiles
        at_mu(now_mu() + 10000)
        # normal signal
        self.res00 = self.urukul1_ch1.set_mu(self.freq_ftw0, asf=self.ampl_asf0, profile=0,
                                             pow_=self.phase_ch0_offset_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # todo: further adjustment, ensure change by sine
        self.res01 = self.urukul1_ch1.set_mu(self.freq_ftw0, asf=self.ampl_asf0, profile=1,
                                             pow_=self.phase_ch0_offset_pow, phase_mode=PHASE_MODE_CONTINUOUS)

        # reference signal
        self.res10 = self.urukul1_ch2.set_mu(self.freq_ftw1, asf=self.ampl_asf0, profile=0,
                                             pow_=self.phase_ch1_offset_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        self.res11 = self.urukul1_ch2.set_mu(self.freq_ftw1, asf=self.ampl_asf0, profile=0,
                                             pow_=self.phase_ch1_offset_pow, phase_mode=PHASE_MODE_CONTINUOUS)

        # set up registers
        at_mu(now_mu() + 10000)
        with parallel:
            self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16))
            self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16))
        with parallel:
            self.urukul1_ch1.set_cfr2(matched_latency_enable=1)
            self.urukul1_ch2.set_cfr2(matched_latency_enable=1)

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



        # cleanup
        with parallel:
            self.urukul1_cpld.cfg_switches(0b0000)
            self.ttl8.off()
            self.ttl9.off()
        self.core.reset()

    def analyze(self):
        print('\n\tresults:')
        print('\t\turukul0 profile 0 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res00)))
        print('\t\turukul1 profile 0 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res10)))
        print('\n')
        print('\t\turukul0 profile 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res01)))
        print('\t\turukul1 profile 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.res11)))
        print('\n')
