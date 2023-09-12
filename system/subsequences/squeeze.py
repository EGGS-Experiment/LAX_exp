from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Squeeze(LAXSubsequence):
    """
    Subsequence: Squeeze

    Squeeze and unsqueeze the ion by applying a quadrupole tone at 2x the secular frequency.
    """
    name = 'squeeze'


    def build_subsequence(self):
        self.setattr_argument('freq_squeeze_khz',           NumberValue(default=10000, ndecimals=3, step=10, min=1, max=400000), group='squeeze')
        self.setattr_argument('time_squeeze_us',            NumberValue(default=20., ndecimals=3, step=100, min=1, max=1000000), group='squeeze')
        self.setattr_argument('att_squeeze_db',             NumberValue(default=10., ndecimals=1, step=0.5, min=0, max=31.5), group='squeeze')
        self.setattr_argument('phase_squeeze_turns',        NumberValue(default=0., ndecimals=3, step=0.1, min=-1., max=1.), group='squeeze')
        self.setattr_argument('phase_antisqueeze_turns',    NumberValue(default=0., ndecimals=3, step=0.1, min=-1., max=1.), group='squeeze')

        # get relevant devices
        self.setattr_device('dds_modulation')
        # tmp remove
        self.setattr_device('ttl8')
        # tmp remove

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.freq_squeeze_ftw =                             self.dds_modulation.frequency_to_ftw(self.freq_squeeze_khz * kHz)
        self.time_squeeze_mu =                              self.core.seconds_to_mu(self.time_squeeze_us * us)
        self.att_squeeze_mu =                               att_to_mu(self.att_squeeze_db * dB)
        self.phase_squeeze_pow =                            self.dds_modulation.turns_to_pow(self.phase_squeeze_turns + 0.)
        self.phase_antisqueeze_pow =                        self.dds_modulation.turns_to_pow(self.phase_antisqueeze_turns + 0.)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        self.core.break_realtime()

        # set dds attenuation here - ensures that dds channel will have correct attenuation
        # for any sequences recorded into DMA during initialize_experiment
        self.dds_modulation.set_att_mu(self.att_squeeze_mu)

        # set up DDS to reinitialize phase each time we set waveform values
        self.dds_modulation.set_mu(self.freq_squeeze_ftw, asf=self.dds_modulation.ampl_modulation_asf, pow_=self.phase_squeeze_pow, profile=0)
        self.dds_modulation.set_mu(self.freq_squeeze_ftw, asf=self.dds_modulation.ampl_modulation_asf, pow_=self.phase_antisqueeze_pow, profile=1)
        self.dds_modulation.reset_phase()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def squeeze(self):
        # reset DDS phase and wait for reset to latch
        self.dds_modulation.set_profile(0)
        at_mu(now_mu() & ~0x7)
        self.dds_modulation.io_update()

        # squeeze for given time
        self.ttl8.on()
        self.dds_modulation.on()
        delay_mu(self.time_squeeze_mu)
        self.dds_modulation.off()
        self.ttl8.off()


    @kernel(flags={"fast-math"})
    def antisqueeze(self):
        # reset DDS phase and wait for reset to latch
        # with parallel:
        # self.dds_modulation.on()
        self.dds_modulation.set_profile(1)
        # with sequential:
        # self.dds_modulation.reset_phase()
        # delay_mu(TIME_PHASEAUTOCLEAR_DELAY_MU)

        # antisqueeze for given time
        self.ttl8.on()
        self.dds_modulation.on()
        delay_mu(self.time_squeeze_mu)
        self.dds_modulation.off()
        self.ttl8.off()

    @kernel(flags={"fast-math"})
    def run(self):
        pass
