from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class PhaserPulse(LAXSubsequence):
    """
    Subsequence: Phaser Pulse

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    """
    name = 'phaser_pulse'

    def build_subsequence(self):
        self.setattr_device('phaser0')

    def prepare_subsequence(self):
        pass
        # self.time_rabiflop_mu = self.get_parameter('time_rabiflop_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        pass

    @kernel(flags={"fast-math"})
    def run(self):
        pass
        # # population transfer pulse
        # self.qubit.on()
        # delay_mu(self.time_rabiflop_mu)
        # self.qubit.off()
