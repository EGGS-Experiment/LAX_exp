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
        self.setattr_argument("enable_squeezing",           BooleanValue(default=True), group=self.name)
        self.setattr_argument("enable_antisqueezing",       BooleanValue(default=True), group=self.name)
        self.setattr_argument('att_squeeze_db',             NumberValue(default=13., ndecimals=1, step=0.5, min=0, max=31.5), group=self.name)

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
        self.time_squeeze_mu =                              np.int64(0)

        # configure squeezing & antisqueezing on/off
        self.squeeze_func =                                 self.dds_modulation.on
        self.antisqueeze_func =                             self.dds_modulation.on
        if not self.enable_squeezing:
            self.squeeze_func =                             self.dds_modulation.off
        if  not self.enable_antisqueezing:
            self.antisqueeze_func =                         self.dds_modulation.off

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set dds attenuation here - ensures that dds channel will have correct attenuation
        # for any sequences recorded into DMA during initialize_experiment
        self.dds_modulation.set_att_mu(self.att_squeeze_mu)
        # self.urukul1_ch2.set_att_mu(self.att_squeeze_mu)
        self.urukul1_ch2.set_att_mu(0x00)


    @kernel(flags={"fast-math"})
    def squeeze(self):
        # set blank output waveform
        self.dds_modulation.set_profile(2)

        # ensure phase_autoclear is enabled ahead of time
        self.urukul1_ch2.write32(_AD9910_REG_CFR1,
                                 ((1 << 16) |
                                  (1 << 13)))
        self.dds_modulation.write32(_AD9910_REG_CFR1,
                                    (1 << 16) | # select_sine_output
                                    (1 << 13))  # phase_autoclear

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7
        # begin output waveform
        at_mu(time_start_mu)
        self.dds_modulation.set_profile(0)

        # open RF switches early since they have a ~100 ns rise time
        at_mu(time_start_mu
              + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              - TIME_URUKUL_RFSWITCH_DELAY_MU)
        self.squeeze_func()
        # self.urukul1_ch2.sw.on()
        self.urukul1_ch2.sw.off()

        # send debug trigger when waveform begins
        # note: 20 is a fudge factor to compensate for stuff (probably things like pipeline delays?)
        # all i know is the fudge is 20ns for ttl to synchronize with DDS signals on oscope
        at_mu(time_start_mu
              + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              + 20)
        self.ttl8.on()

        # squeeze for given time
        # note: turn off switch early to account for switch rise time
        delay_mu(self.time_squeeze_mu
                 - TIME_ZASWA2_SWITCH_DELAY_MU)
        self.dds_modulation.off()
        # self.urukul1_ch2.sw.on()
        self.urukul1_ch2.sw.on()

        # send debug trigger after switch has fully closed
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)
        self.ttl8.off()

        # unset phase_autoclear to keep output coherent
        self.dds_modulation.write32(_AD9910_REG_CFR1, (1 << 16))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16))


    @kernel(flags={"fast-math"})
    def antisqueeze(self):
        # unset phase_autoclear to keep output coherent
        self.dds_modulation.write32(_AD9910_REG_CFR1, (1 << 16))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16))

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
        self.antisqueeze_func()
        # self.urukul1_ch2.sw.on()
        self.urukul1_ch2.sw.off()

        # send debug trigger when waveform begins
        at_mu(time_start_mu
              + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              + 20)
        self.ttl9.on()

        # squeeze for given time
        # note: turn off switch early to account for switch rise time
        delay_mu(self.time_squeeze_mu - TIME_ZASWA2_SWITCH_DELAY_MU)
        self.dds_modulation.off()
        self.urukul1_ch2.sw.off()

        # send debug trigger after switch has fully closed
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)
        self.ttl9.off()

    @kernel(flags={"fast-math"})
    def configure(self, freq_ftw: TInt32, phase_pow: TInt32, time_mu: TInt64) -> TInt64:
        # calculate squeezing frequency half-period
        time_half_period_ns =                               1. / (2. * ns * self.dds_modulation.ftw_to_frequency(freq_ftw))

        # ensure squeezing period is close to a multiple of the frequency half-period
        self.time_squeeze_mu =                              time_mu
        if self.time_squeeze_mu % time_half_period_ns:
            # round squeezing time up to the nearest multiple of the frequency half-period
            num_half_period_multiples =                     round(self.time_squeeze_mu / time_half_period_ns + 0.5)
            self.time_squeeze_mu =                          np.int64(round(time_half_period_ns * num_half_period_multiples))

        # calculate antisqueeze phase value
        # phase_antisqueeze_turns =                           (0.5
        #                                                      - ((self.time_squeeze_mu / (2. * time_half_period_ns)) % 1)
        #                                                      - self.dds_modulation.pow_to_turns(phase_pow))
        phase_antisqueeze_turns =                           0.5 + self.dds_modulation.pow_to_turns(phase_pow)
        phase_antisqueeze_pow =                             self.dds_modulation.turns_to_pow(phase_antisqueeze_turns)
        # print(time_half_period_ns)
        # print(self.time_squeeze_mu)
        # print(phase_antisqueeze_turns)
        # self.core.break_realtime()

        # set waveforms for profiles
        at_mu(now_mu() + 10000)
        # set up DDS to track phase when we change profiles
        # squeezing waveform
        self.dds_modulation.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        # antisqueezing waveform
        self.dds_modulation.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=1,
                                   pow_=phase_antisqueeze_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # blank waveform
        self.dds_modulation.set_mu(freq_ftw, asf=0x0, profile=2,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)

        # at_mu(now_mu() + 10000)
        # self.urukul1_ch2.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0,
        #                            pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        # self.urukul1_ch2.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=1,
        #                            pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        # self.urukul1_ch2.set_mu(freq_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=2,
        #                            pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(freq_ftw, asf=0x0, profile=0,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(freq_ftw, asf=0x0, profile=1,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(freq_ftw, asf=0x0, profile=2,
                                   pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)

        # return time value
        return self.time_squeeze_mu

    @kernel(flags={"fast-math"})
    def run(self):
        pass
