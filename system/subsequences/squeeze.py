from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Squeeze(LAXSubsequence):
    """
    Subsequence: Squeeze

    Squeeze and unsqueeze the ion by applying a quadrupole tone at 2x the secular frequency.
    """
    name = 'squeeze'


    def build_subsequence(self):
        self.setattr_argument('freq_squeeze_khz',           NumberValue(default=2176, ndecimals=3, step=10, min=1, max=400000), group='squeeze')
        self.setattr_argument('att_squeeze_db',             NumberValue(default=10., ndecimals=1, step=0.5, min=0, max=31.5), group='squeeze')

        self.setattr_argument('time_squeeze_us',            NumberValue(default=5., ndecimals=3, step=100, min=1, max=1000000), group='squeeze')
        self.setattr_argument('phase_squeeze_turns',        NumberValue(default=0., ndecimals=3, step=0.1, min=-1., max=1.), group='squeeze')
        self.setattr_argument('phase_antisqueeze_turns',    NumberValue(default=0., ndecimals=3, step=0.1, min=-1., max=1.), group='squeeze')

        # get relevant devices
        self.setattr_device('dds_modulation')
        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        # tmp remove

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.freq_squeeze_ftw =                             self.dds_modulation.frequency_to_ftw(self.freq_squeeze_khz * kHz)
        self.att_squeeze_mu =                               att_to_mu(self.att_squeeze_db * dB)

        self.time_squeeze_mu =                              self.core.seconds_to_mu(self.time_squeeze_us * us)
        self.phase_squeeze_pow =                            self.dds_modulation.turns_to_pow(self.phase_squeeze_turns + 0.)
        self.phase_antisqueeze_pow =                        self.dds_modulation.turns_to_pow(self.phase_antisqueeze_turns + 0.5)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        self.core.break_realtime()

        # set dds attenuation here - ensures that dds channel will have correct attenuation
        # for any sequences recorded into DMA during initialize_experiment
        self.dds_modulation.set_att_mu(self.att_squeeze_mu)

        # set up DDS to track phase when we change profiles
        # squeezing waveform
        self.dds_modulation.set_mu(self.freq_squeeze_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0,
                                   pow_=self.phase_squeeze_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # antisqueezing waveform
        self.dds_modulation.set_mu(self.freq_squeeze_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=1,
                                   pow_=self.phase_antisqueeze_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # blank waveform
        self.dds_modulation.set_mu(self.freq_squeeze_ftw, asf=0x0, profile=2,
                                   pow_=self.phase_antisqueeze_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def squeeze(self):
        # set blank output waveform
        self.dds_modulation.set_profile(2)
        # ensure phase_autoclear is enabled ahead of time
        self.dds_modulation.write32(_AD9910_REG_CFR1,
                                    (1 << 16) | # select_sine_output
                                    (1 << 13))  # phase_autoclear

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7

        # begin output waveform
        # note: no need to reset phase since phase_autoclear flag should be set
        # in initialize_subsequence, causing phase accumulator reset upon profile change
        at_mu(time_start_mu)
        self.dds_modulation.set_profile(0)

        # open RF switches early since they have a ~100 ns rise time
        at_mu(time_start_mu +
              (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              - TIME_URUKUL_RFSWITCH_DELAY_MU)
        self.dds_modulation.on()

        # tmp remove
        # send debug trigger when waveform begins
        at_mu(time_start_mu + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU))
        self.ttl8.on()
        # tmp remove

        # squeeze for given time
        delay_mu(self.time_squeeze_mu)
        self.ttl8.off()
        self.dds_modulation.off()


    @kernel(flags={"fast-math"})
    def antisqueeze(self):
        # unset phase_autoclear flag to ensure phase remains tracked,
        # and ensure set_sine_output flag remains set
        self.dds_modulation.write32(_AD9910_REG_CFR1, (1 << 16))

        # align to coarse RTIO clock to prepare output waveform for antisqueezing
        at_mu(now_mu() & ~0x7)
        self.dds_modulation.set_profile(1)

        # enable antisqueezing output for given time
        self.ttl9.on()
        self.dds_modulation.on()
        delay_mu(self.time_squeeze_mu)

        # stop antisqueezing output
        self.ttl9.off()
        self.dds_modulation.off()

    @kernel(flags={"fast-math"})
    def run(self):
        pass
