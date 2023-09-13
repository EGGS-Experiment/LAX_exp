import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SqueezeConfigurable(LAXSubsequence):
    """
    Subsequence: Squeeze (Configurable)

    Squeeze and unsqueeze the ion by applying a quadrupole tone at 2x the secular frequency.
    Allows configuration of the squeezing to sweep over the values.
    """
    name = 'squeeze_configurable'


    def build_subsequence(self):
        self.setattr_argument('att_squeeze_db',             NumberValue(default=10., ndecimals=1, step=0.5, min=0, max=31.5), group='squeeze_configurable')
        self.setattr_argument("enable_antisqueezing",       BooleanValue(default=True), group='squeeze_configurable')

        # get relevant devices
        self.setattr_device('dds_modulation')
        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('urukul1_ch2')
        # tmp remove

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.att_squeeze_mu =                               att_to_mu(self.att_squeeze_db * dB)

        # create empty holder variables to support later configuration
        self.freq_squeeze_ftw =                             np.int32(0)
        self.phase_squeeze_pow =                            np.int32(0)
        self.phase_antisqueeze_pow =                        np.int32(0)

        # configure antisqueezing on/off
        if self.enable_antisqueezing:
            self.antisqueeze_func =                         self.dds_modulation.on
        else:
            self.antisqueeze_func =                         self.dds_modulation.off

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set dds attenuation here - ensures that dds channel will have correct attenuation
        # for any sequences recorded into DMA during initialize_experiment
        self.dds_modulation.set_att_mu(self.att_squeeze_mu)

        # tmp remove
        self.urukul1_ch2.set_att_mu(self.att_squeeze_mu)
        # tmp remove


    @kernel(flags={"fast-math"})
    def squeeze(self):
        # set blank output waveform
        self.dds_modulation.set_profile(2)

        # ensure phase_autoclear is enabled ahead of time
        # tmp remove
        self.urukul1_ch2.write32(_AD9910_REG_CFR1,
                                 ((1 << 16) |
                                  (1 << 13)))
        # tmp remove
        self.dds_modulation.write32(_AD9910_REG_CFR1,
                                    (1 << 16) | # select_sine_output
                                    (1 << 13))  # phase_autoclear

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7
        # begin output waveform
        at_mu(time_start_mu)
        self.dds_modulation.set_profile(0)

        # open RF switches early since they have a ~100 ns rise time
        at_mu(time_start_mu +
              (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              - TIME_URUKUL_RFSWITCH_DELAY_MU)
        self.dds_modulation.on()
        # tmp remove
        self.urukul1_ch2.sw.on()
        # tmp remove

        # tmp remove
        # send debug trigger when waveform begins
        at_mu(time_start_mu + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU))
        self.ttl8.on()
        # tmp remove

        # # squeeze for given time
        # delay_mu(self.time_squeeze_mu)
        # self.dds_modulation.off()
        # # tmp remove
        # self.ttl8.off()
        # self.urukul1_ch2.sw.off()
        # # tmp remove


    @kernel(flags={"fast-math"})
    def antisqueeze(self):
        # unset phase_autoclear flag to ensure phase remains tracked,
        # and ensure set_sine_output flag remains set
        # tmp remove
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        # tmp remove
        self.dds_modulation.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        # set blank profile to ensure switching is exact
        self.dds_modulation.set_profile(2)

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7
        # begin output waveform
        at_mu(time_start_mu)
        self.dds_modulation.set_profile(1)

        # open RF switches early since they have a ~100 ns rise time
        at_mu(time_start_mu +
              (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              - TIME_URUKUL_RFSWITCH_DELAY_MU)
        # self.dds_modulation.on()
        self.antisqueeze_func()
        # tmp remove
        self.urukul1_ch2.sw.on()
        # tmp remove

        # tmp remove
        # send debug trigger when waveform begins
        at_mu(time_start_mu + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU))
        self.ttl9.on()
        # tmp remove

        # # antisqueeze for given time
        # delay_mu(self.time_squeeze_mu)
        # self.dds_modulation.off()
        # # tmp remove
        # self.ttl9.off()
        # self.urukul1_ch2.sw.off()
        # # tmp remove

    @kernel(flags={"fast-math"})
    def configure(self, freq_ftw: TInt32, phase_pow: TInt32):
        # calculate antisqueeze phase value
        phase_antisqueeze_pow = self.dds_modulation.turns_to_pow(self.dds_modulation.pow_to_turns(phase_pow)
                                                                 + 0.5)

        # set waveforms for profiles
        at_mu(now_mu() + 5000)
        # set up DDS to track phase when we change profiles
        # squeezing waveform
        self.dds_modulation.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0,
                                   pow_=phase_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # antisqueezing waveform
        self.dds_modulation.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=1,
                                   pow_=phase_antisqueeze_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # blank waveform
        self.dds_modulation.set_mu(freq_ftw, asf=0x0, profile=2,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.core.break_realtime()

        # tmp remove
        self.urukul1_ch2.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=1,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(freq_ftw, asf=0x0, profile=2,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)
        # tmp remove

    @kernel(flags={"fast-math"})
    def run(self):
        pass
