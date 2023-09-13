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
        self.frequency_mhz =                100
        self.amplitude_pct =                50.
        self.att_db =                       10.

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
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        # prepare device parameter values
        self._prepare_parameters()

    def _prepare_parameters(self):
        # calculate base waveform values
        self.ampl_asf0 =                    self.urukul1_ch1.amplitude_to_asf(self.amplitude_pct / 100.)
        self.ampl_asf1 =                    self.urukul1_ch1.amplitude_to_asf(self.amplitude_pct / 100. * 1.5)

        self.freq_ftw0 =                    self.urukul1_ch1.frequency_to_ftw(self.frequency_mhz * MHz)
        self.freq_ftw1 =                    self.urukul1_ch1.frequency_to_ftw(self.frequency_mhz * MHz)

        self.phase_ch0_offset_pow =         self.urukul1_ch1.turns_to_pow(self.phase_turns0)
        self.phase_ch1_offset_pow =         self.urukul1_ch1.turns_to_pow(self.phase_turns1)

        # calculate phase and timing values to compensate for ch0 & ch1 delays
        self.phase_ch1_inherent_turns =     0.472

        # self.time_ch1_latency_ns =          12.13
        self.time_ch1_latency_ns =          0.
        self.phase_ch1_latency_turns =      (self.frequency_mhz * MHz) * (self.time_ch1_latency_ns * ns)
        self.time_delay_mu =                self.core.seconds_to_mu(self.time_delay_ns * ns)

        # combine compensation values into ch1 phase: ch0 vs ch1 inherent phase shift/relation, ch0 vs ch1 inherent time delay, and ch1 offset
        self.phase_ch1_final_pow =          self.urukul1_ch1.turns_to_pow(self.phase_ch1_inherent_turns +
                                                                          self.phase_ch1_latency_turns +
                                                                          self.phase_turns1)
        self.phase_pi_pow =                 self.urukul1_ch1.turns_to_pow(0.5)

        # tmp remove
        self.res00 = np.int32(0)
        self.res10 = np.int32(0)

        self.res01 = np.int32(0)
        self.res11 = np.int32(0)

        self.t00 = np.int64(0)
        self.t01 = np.int64(0)
        # tmp remove

    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        # prepare asynchronous/timing-agnostic stuff
        # with parallel:
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
                                             pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.res01 = self.urukul1_ch1.set_mu(self.freq_ftw0, asf=self.ampl_asf0, profile=1,
                                             pow_=self.phase_pi_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch1.set_mu(self.freq_ftw0, asf=0x0, profile=2,
                                pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)

        # reference signal
        self.res10 = self.urukul1_ch2.set_mu(self.freq_ftw1, asf=self.ampl_asf0, profile=0,
                                             pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.res11 = self.urukul1_ch2.set_mu(self.freq_ftw1, asf=self.ampl_asf0, profile=1,
                                             pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(self.freq_ftw1, asf=0x0, profile=2,
                                pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)

        # set up registers
        at_mu(now_mu() + 10000)
        self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)

        # turn switches off and set a 0 amplitude waveform
        at_mu(now_mu() + 10000)
        with parallel:
            self.urukul1_cpld.set_profile(2)
            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        '''PREPARE DEVICES AND TIMING'''
        self._run_prepare()

        delay_mu(125000)
        # note: we coarse align to previous SYNC_CLK period
        time_start_mu = now_mu() & ~7


        '''SQUEEZE START'''
        # set start (active) profile
        at_mu(time_start_mu)
        self.urukul1_cpld.set_profile(0)

        # open RF switches early since they have ~100 ns rise time
        at_mu(time_start_mu + ((416 + 63) - 140))
        self.urukul1_ch1.sw.on()
        self.urukul1_ch2.sw.on()

        # send trigger when waveform begins
        at_mu(time_start_mu + (416 + 63))
        self.ttl8.on()
        self.t00 = now_mu()


        '''SQUEEZE STOP'''
        delay_mu(self.time_delay_mu)
        self.t01 = now_mu()

        self.ttl8.off()
        self.urukul1_ch1.sw.off()
        # self.urukul1_ch2.sw.off()


        '''
        ARBITRARY PULSE SEQUENCE
        '''
        # self.t00 = now_mu()
        with parallel:
            delay_mu(10000)
            with sequential:
                self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16))
                self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16))
                self.urukul1_cpld.set_profile(1)
        # self.t01 = now_mu()


        '''ANTISQUEEZE'''
        time_antisqueeze_mu = now_mu()

        # open RF switch in advance
        at_mu(time_antisqueeze_mu - 140)
        self.urukul1_ch1.sw.on()

        at_mu(time_antisqueeze_mu)
        self.ttl9.on()
        delay_mu(1000)
        self.ttl9.off()


        '''CLEANUP'''
        with parallel:
            self.ttl8.off()
            self.ttl9.off()
            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()
        self.core.reset()

    def analyze(self):
        print('\n\tresults:')
        print('\t\turukul1_ch1 profile 0 phase: {:.4f}'.format(self.urukul1_ch1.pow_to_turns(self.res00)))
        print('\t\turukul1_ch2 profile 0 phase: {:.4f}'.format(self.urukul1_ch1.pow_to_turns(self.res10)))
        print('\n')
        print('\t\turukul1_ch1 profile 1 phase: {:.4f}'.format(self.urukul1_ch1.pow_to_turns(self.res01)))
        print('\t\turukul1_ch2 profile 1 phase: {:.4f}'.format(self.urukul1_ch1.pow_to_turns(self.res11)))
        print('\n')
        print('\t\tt00: {:d}'.format(self.t00))
        print('\t\tt01: {:d}'.format(self.t01))
        print('\t\tdelay: {:d}'.format(self.t01 - self.t00))

