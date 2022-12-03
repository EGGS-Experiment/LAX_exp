from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class Tickle(LAXSubsequence):
    """
    Subsequence: Tickle
        Apply the tickle beam for a given time.
    """
    name = 'tickle'

    devices = [
        'tickle'
    ]
    subsequence_parameters = {
        'time_tickle_mu':           ('timing.time_tickle_mu', us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # readout pulse
        self.tickle.cfg_sw(1)
        delay_mu(self.time_tickle_mu)
        self.tickle.cfg_sw(0)
