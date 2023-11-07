import numpy as np
from artiq.experiment import *

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1


class UrukulLatencyTest(EnvExperiment):
    """
    Urukul Latency Test
    Test the relative output latency between different Urukul boards and channels.
    """

    def build(self):
        # waveform
        self.freq_mhz =                     151.2
        self.ampl_pct =                     50.

        # attenuation
        self.att_db =                           6.

        # time
        self.time_pulse_us =                    10000.

        # phase - global
        self.phase_rsb_turns =                  0.
        self.phase_bsb_turns =                  0.
        self.phase_carrier_turns =              0.

        self.phase_urukul1_adjust_turns =       0.
        self.time_urukul1_system_latency_ns =   0.

        # phase - channel
        self.time_urukul0_ch2_latency_ns =      1.91
        self.time_urukul0_ch3_latency_ns =      1.97

        # self.time_urukul1_ch1_latency_ns =      -0.44
        # self.time_urukul1_ch2_latency_ns =      -1.41
        # self.time_urukul1_ch3_latency_ns =      -0.34
        self.time_urukul1_ch1_latency_ns =      0.
        self.time_urukul1_ch2_latency_ns =      0.
        self.time_urukul1_ch3_latency_ns =      0.

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
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")
        self.setattr_device("urukul2_ch2")
        self.setattr_device("urukul2_ch3")

        self.dds_cpld = self.urukul0_cpld
        self.dds = self.urukul0_ch1

        # prepare device parameter values
        self._prepare_parameters_general()
        self._prepare_parameters_phase()

    def _prepare_parameters_general(self):
        # attenuation
        self.att_mu =                       self.dds_cpld.att_to_mu(self.att_db * dB)
        # preallocate register storage for urukul
        self._reg_urukul0_current =         np.int32(0)
        self._reg_urukul1_current =         np.int32(0)

        # amplitude
        self.ampl_rsb_asf =                 self.dds.amplitude_to_asf(self.ampl_rsb_pct / 100.)
        self.ampl_bsb_asf =                 self.dds.amplitude_to_asf(self.ampl_bsb_pct / 100.)
        self.ampl_carrier_asf =             self.dds.amplitude_to_asf(self.ampl_carrier_pct / 100.)

        # pulse time
        self.time_pulse_mu =                self.core.seconds_to_mu(self.time_pulse_us * us)

        # empty waveform
        self.freq_empty_ftw =               self.dds.frequency_to_ftw(228.12 * MHz)

    def _prepare_parameters_phase(self):
        # calculate frequency values
        freq_rsb_hz =                       self.freq_carrier_mhz * MHz - self.freq_sideband_khz * kHz + self.freq_offset_khz * kHz
        freq_bsb_hz =                       self.freq_carrier_mhz * MHz + self.freq_sideband_khz * kHz + self.freq_offset_khz * kHz
        freq_carrier_hz =                   self.freq_carrier_mhz * MHz
        # freq_rsb_hz =                       self.freq_carrier_mhz * MHz
        # freq_bsb_hz =                       self.freq_carrier_mhz * MHz
        # freq_carrier_hz =                   self.freq_carrier_mhz * MHz

        self.freq_rsb_ftw =                 self.dds.frequency_to_ftw(freq_rsb_hz)
        self.freq_bsb_ftw =                 self.dds.frequency_to_ftw(freq_bsb_hz)
        self.freq_carrier_ftw =             self.dds.frequency_to_ftw(freq_carrier_hz)

        # calculate phase values for each channel
        self.phase_urukul0_ch1_turns =      self.phase_rsb_turns
        self.phase_urukul0_ch2_turns =      freq_bsb_hz * (self.time_urukul0_ch2_latency_ns * ns) + self.phase_bsb_turns
        self.phase_urukul0_ch3_turns =      freq_carrier_hz * (self.time_urukul0_ch3_latency_ns * ns) + self.phase_carrier_turns

        self.phase_urukul1_ch1_turns =      (freq_rsb_hz * ((self.time_urukul1_ch1_latency_ns - self.time_urukul1_system_latency_ns) * ns)
                                             + self.phase_urukul1_adjust_turns
                                             + self.phase_rsb_turns)
        self.phase_urukul1_ch2_turns =      (freq_bsb_hz * ((self.time_urukul1_ch2_latency_ns - self.time_urukul1_system_latency_ns) * ns)
                                             + self.phase_urukul1_adjust_turns
                                             + self.phase_bsb_turns)
        self.phase_urukul1_ch3_turns =      (freq_carrier_hz * ((self.time_urukul1_ch3_latency_ns - self.time_urukul1_system_latency_ns) * ns)
                                             + self.phase_urukul1_adjust_turns
                                             + self.phase_carrier_turns
                                             + 0.)

        self.phase_urukul0_ch1_pow =        self.dds.turns_to_pow(self.phase_urukul0_ch1_turns)
        self.phase_urukul0_ch2_pow =        self.dds.turns_to_pow(self.phase_urukul0_ch2_turns)
        self.phase_urukul0_ch3_pow =        self.dds.turns_to_pow(self.phase_urukul0_ch3_turns)

        self.phase_urukul1_ch1_pow =        self.dds.turns_to_pow(self.phase_urukul1_ch1_turns)
        self.phase_urukul1_ch2_pow =        self.dds.turns_to_pow(self.phase_urukul1_ch2_turns)
        self.phase_urukul1_ch3_pow =        self.dds.turns_to_pow(self.phase_urukul1_ch3_turns)


    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        """
        Prepare asynchronous/timing-agnostic stuff.
        """
        '''PREPARE - TTLs'''
        self.ttl8.off()
        self.ttl9.off()
        delay_mu(100)
        with parallel:
            self.urukul0_ch1.sw.off()
            self.urukul0_ch2.sw.off()
            self.urukul0_ch3.sw.off()

            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()
            self.urukul1_ch3.sw.off()


        '''PREPARE - DDSs'''
        at_mu(now_mu() + 25000)
        self._reg_urukul0_current = self.urukul0_cpld.get_att_mu()
        self._reg_urukul0_current &= (0xFF << 0)
        self._reg_urukul0_current |= ((self.att_mu << 8) | (self.att_mu << 16) | (self.att_mu << 24))
        at_mu(now_mu() + 10000)
        # self.urukul0_cpld.set_all_att_mu(self._reg_urukul0_current)
        self.urukul0_ch1.set_att_mu(self.att_mu)
        at_mu(now_mu() + 10000)
        self.urukul0_ch2.set_att_mu(self.att_mu)
        at_mu(now_mu() + 10000)
        self.urukul0_ch3.set_att_mu(self.att_mu)


        at_mu(now_mu() + 25000)
        self._reg_urukul1_current = self.urukul1_cpld.get_att_mu()
        self._reg_urukul1_current &= (0xFF << 0)
        self._reg_urukul1_current |= ((self.att_mu << 8) | (self.att_mu << 16) | (self.att_mu << 24))
        at_mu(now_mu() + 10000)
        self.urukul1_cpld.set_all_att_mu(self._reg_urukul1_current)

        at_mu(now_mu() + 25000)
        self.urukul0_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)

        self.urukul1_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS WAVEFORMS - ACTIVE'''
        at_mu(now_mu() + 25000)
        self.urukul0_ch1.set_mu(self.freq_rsb_ftw, asf=self.ampl_rsb_asf, pow_=self.phase_urukul0_ch1_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_mu(self.freq_bsb_ftw, asf=self.ampl_bsb_asf, pow_=self.phase_urukul0_ch2_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, pow_=self.phase_urukul0_ch3_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)

        at_mu(now_mu() + 25000)
        self.urukul1_ch1.set_mu(self.freq_rsb_ftw, asf=self.ampl_rsb_asf, pow_=self.phase_urukul1_ch1_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(self.freq_bsb_ftw, asf=self.ampl_bsb_asf, pow_=self.phase_urukul1_ch2_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, pow_=self.phase_urukul1_ch3_pow,
                                profile=0, phase_mode=PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS WAVEFORMS - EMPTY'''
        at_mu(now_mu() + 25000)
        self.urukul0_ch1.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)

        at_mu(now_mu() + 25000)
        self.urukul1_ch1.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_mu(self.freq_empty_ftw, asf=0x01, pow_=0x01, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS REGISTERS'''
        at_mu(now_mu() + 25000)
        self.urukul0_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 25000)
        self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 25000)
        self.urukul0_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch3.set_cfr2(matched_latency_enable=1)

        at_mu(now_mu() + 25000)
        self.urukul1_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch3.set_cfr2(matched_latency_enable=1)


        '''PREPARE - INITIALIZATION'''
        at_mu(now_mu() + 10000)
        with parallel:
            self.urukul0_cpld.set_profile(1)
            self.urukul1_cpld.set_profile(1)

        at_mu(now_mu() + 1000)
        with parallel:
            self.urukul0_ch1.sw.on()
            self.urukul0_ch2.sw.on()
            self.urukul0_ch3.sw.on()

            self.urukul1_ch1.sw.on()
            self.urukul1_ch2.sw.on()
            self.urukul1_ch3.sw.on()


    @kernel(flags={"fast-math"})
    def run(self):
        """
        todo: document
        """
        '''PREPARE'''
        self.core.reset()
        self._run_prepare()

        delay_mu(125000)
        # note: we coarse align to previous SYNC_CLK period
        time_start_mu = now_mu() & ~0x7


        '''PULSE START'''
        # set start (active) profile
        at_mu(time_start_mu)
        self.urukul0_cpld.set_profile(0)

        # open RF switches early since they have ~100 ns rise time
        # at_mu(time_start_mu + ((416 + 63) - 200))
        # with parallel:
        #     self.urukul0_ch1.sw.on()
        #     self.urukul0_ch2.sw.on()
        #     self.urukul0_ch3.sw.on()
        #
        #     self.urukul1_ch1.sw.on()
        #     self.urukul1_ch2.sw.on()
        #     self.urukul1_ch3.sw.on()

        # send trigger when waveform begins
        at_mu(time_start_mu + (416 + 63))
        self.ttl8.on()
        delay_mu(self.time_pulse_mu)


        '''PULSE STOP'''
        # close RF switches
        with parallel:
            self.urukul0_ch1.sw.off()
            self.urukul0_ch2.sw.off()
            self.urukul0_ch3.sw.off()

            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()
            self.urukul1_ch3.sw.off()

        self.ttl8.off()


    def analyze(self):
        print('\n\tresults:')
        print('\t\turukul0 channel 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul0_ch1_pow)))
        print('\t\turukul1 channel 1 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul1_ch1_pow)))
        print('\n')
        print('\t\turukul0 channel 2 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul0_ch2_pow)))
        print('\t\turukul1 channel 2 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul1_ch2_pow)))
        print('\n')
        print('\t\turukul0 channel 3 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul0_ch3_pow)))
        print('\t\turukul1 channel 3 phase: {:.4f}'.format(self.urukul0_ch0.pow_to_turns(self.phase_urukul1_ch3_pow)))
        print('\n')
        print(self.dds.ftw_to_frequency(self.freq_rsb_ftw) / MHz)
        print(self.dds.ftw_to_frequency(self.freq_bsb_ftw) / MHz)
        print(self.dds.ftw_to_frequency(self.freq_carrier_ftw) / MHz)
