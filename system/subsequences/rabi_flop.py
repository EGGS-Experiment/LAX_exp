from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
# todo: make rabi flopping time an argument instead of parameter and require it to be specified each time


class RabiFlop(LAXSubsequence):
    """
    Subsequence: Rabi Flop

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    """
    name = 'rabi_flop'
    kernel_invariants = {
        "time_rabiflop_mu"
    }

    def build_subsequence(self):
        self.setattr_device('qubit')

    def prepare_subsequence(self):
        self.time_rabiflop_mu = self.get_parameter('time_rabiflop_us', group='timing',
                                                   override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # population transfer pulse
        self.qubit.on()
        delay_mu(self.time_rabiflop_mu)
        self.qubit.off()
