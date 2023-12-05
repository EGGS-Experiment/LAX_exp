from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class FockState(LAXSubsequence):
    """
    Subsequence: Fock State

    Create a pure Fock state with n > 0 by transferring population on the blue sideband
    after cooling to the n=0 ground state.
    """
    name = 'fock_state'

    def build_subsequence(self):
        self.setattr_argument("enable_fock_state",            BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("ampl_fock_pulse_mu",             NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')

        self.time_fock_pulse_mu = self.get_parameter('time_rabiflop_us', group='timing',
                                                   override=True, conversion_function=seconds_to_mu, units=us)
        self.ampl_fock_pulse_mu = self.get_parameter('time_rabiflop_us', group='timing',
                                                   override=True, conversion_function=seconds_to_mu, units=us)

        self.setattr_device('qubit')
        # todo: set arguments

    def prepare_subsequence(self):
        self.time_rabiflop_mu = self.get_parameter('time_rabiflop_us', group='timing',
                                                   override=True, conversion_function=seconds_to_mu, units=us)

        # todo: convert arguments

    @kernel(flags={"fast-math"})
    def run(self):
        # population transfer pulse
        self.qubit.on()
        delay_mu(self.time_rabiflop_mu)
        self.qubit.off()
