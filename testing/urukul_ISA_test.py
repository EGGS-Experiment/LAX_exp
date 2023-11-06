import numpy as np
from artiq.experiment import *

from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1


class UrukulISATest(EnvExperiment):
    """
    Urukul ISA Test
    Testing ISA via Urukuls
    """

    def build(self):
        self.freq_carrier_mhz =             82.
        self.freq_sideband_khz =            891.2
        self.freq_offset_khz =              0.

        self.ampl_rsb_pct =                 50.
        self.ampl_bsb_pct =                 50.
        self.ampl_carrier_pct =             50.

        self.att_db =                       10.

        self.phase_rsb_turns =              0.
        self.phase_bsb_turns =              0.
        # self.phase_carrier_turns =          0.

        self.phase_ch1_inherent_turns =     0.
        self.time_ch1_delay_ns =            100.

        self.time_pulse_us =                100000000.

    def prepare(self):
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # TTLs
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        # urukul devices
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        self.dds_cpld = self.urukul0_cpld
        self.dds = self.urukul0_ch1

        # prepare device parameter values
        self._prepare_parameters()

    def _prepare_parameters(self):
        # calculate base waveform values
        self.att_mu =                       self.dds_cpld.att_to_mu(self.att_db * dB)

        self.ampl_rsb_asf =                 self.dds.amplitude_to_asf(self.ampl_rsb_pct / 100.)
        self.ampl_bsb_asf =                 self.dds.amplitude_to_asf(self.ampl_bsb_pct / 100.)
        self.ampl_carrier_asf =             self.dds.amplitude_to_asf(self.ampl_carrier_pct / 100.)

        self.freq_rsb_ftw =                 self.dds.frequency_to_ftw(self.freq_carrier_mhz * MHz - self.freq_sideband_khz * kHz + self.freq_offset_khz * kHz)
        self.freq_bsb_ftw =                 self.dds.frequency_to_ftw(self.freq_carrier_mhz * MHz + self.freq_sideband_khz * kHz + self.freq_offset_khz * kHz)
        self.freq_carrier_ftw =             self.dds.frequency_to_ftw(self.freq_carrier_mhz * MHz)

        self.phase_rsb_pow_ch0 =            self.dds.turns_to_pow(self.phase_rsb_turns)
        self.phase_bsb_pow_ch0 =            self.dds.turns_to_pow(self.phase_bsb_turns)
        self.phase_carrier_pow_ch0 =        self.dds.turns_to_pow(0.)

        self.phase_rsb_pow_ch1 =            self.dds.turns_to_pow(self.phase_rsb_turns + self.phase_ch1_inherent_turns)
        self.phase_bsb_pow_ch1 =            self.dds.turns_to_pow(self.phase_bsb_turns + self.phase_ch1_inherent_turns)
        self.phase_carrier_pow_ch1 =        self.dds.turns_to_pow(self.phase_ch1_inherent_turns)

        self.time_pulse_mu =                self.core.seconds_to_mu(self.time_pulse_us * us)

        # preallocate register storage for urukul
        self._reg_urukul0_current =         np.int32(0)
        self._reg_urukul1_current =         np.int32(0)

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
        self._reg_urukul0_current &= ~(0xFF << 0)
        self._reg_urukul0_current |= ((self.att_mu << 8) | (self.att_mu << 16) | (self.att_mu << 24))

        at_mu(now_mu() + 25000)
        self._reg_urukul1_current = self.urukul1_cpld.get_att_mu()
        self._reg_urukul1_current &= ~(0xFF << 0)
        self._reg_urukul1_current |= ((self.att_mu << 8) | (self.att_mu << 16) | (self.att_mu << 24))

        at_mu(now_mu() + 25000)
        self.urukul0_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)

        self.urukul1_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS WAVEFORMS'''
        at_mu(now_mu() + 10000)
        self.urukul0_ch1.set_mu(self.freq_rsb_ftw, asf=self.ampl_rsb_asf, profile=0, pow_=self.phase_rsb_pow_ch0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_mu(self.freq_bsb_ftw, asf=self.ampl_bsb_asf, profile=0, pow_=self.phase_bsb_pow_ch0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, profile=0, pow_=self.phase_carrier_pow_ch0, phase_mode=PHASE_MODE_CONTINUOUS)

        self.urukul1_ch1.set_mu(self.freq_rsb_ftw, asf=self.ampl_rsb_asf, profile=0, pow_=self.phase_rsb_pow_ch1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(self.freq_bsb_ftw, asf=self.ampl_bsb_asf, profile=0, pow_=self.phase_bsb_pow_ch1, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, profile=0, pow_=self.phase_carrier_pow_ch1, phase_mode=PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS REGISTERS'''
        at_mu(now_mu() + 10000)
        self.urukul0_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 10000)
        self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 10000)
        self.urukul0_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch3.set_cfr2(matched_latency_enable=1)

        at_mu(now_mu() + 10000)
        self.urukul1_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch3.set_cfr2(matched_latency_enable=1)


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
        time_start_mu = now_mu() & ~7


        '''PULSE START'''
        # set start (active) profile
        at_mu(time_start_mu)
        with parallel:
            self.urukul0_cpld.set_profile(0)
            self.urukul1_cpld.set_profile(0)

        # open RF switches early since they have ~100 ns rise time
        at_mu(time_start_mu + ((416 + 63) - 140))
        with parallel:
            self.urukul0_ch1.sw.on()
            self.urukul0_ch2.sw.on()
            self.urukul0_ch3.sw.on()

            self.urukul1_ch1.sw.on()
            self.urukul1_ch2.sw.on()
            self.urukul1_ch3.sw.on()

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
        pass
