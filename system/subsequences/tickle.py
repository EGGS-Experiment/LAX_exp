from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Tickle(LAXSubsequence):
    """
    Subsequence: Tickle
        Heat the ion by applying a tickle at the secular frequency.
    """
    name = 'tickle'

    parameters = {
        'time_tickle_mu':                   ('timing.time_tickle_us',                   us_to_mu)
    }
    devices = [
        'tickle'
    ]
    
    @kernel(flags={"fast-math"})
    def run(self):
        self.tickle.on()
        delay_mu(self.time_tickle_mu)
        self.tickle.off()
