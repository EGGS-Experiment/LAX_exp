from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class RabiFlop(LAXSubsequence):
    """
    Subsequence: Rabi Flop
        Apply the 729nm beam to cause rabi flopping.
    """
    name = 'rabi_flop'

    devices = [
        'qubit'
    ]
    subsequence_parameters = {
        # 'time_tickle_mu':           ('timing.time_tickle_mu', us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # readout pulse
        self.qubit.cfg_sw(1)
        # delay_mu(self.time_tickle_mu)
        # self.qubit.cfg_sw(0)
