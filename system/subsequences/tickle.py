from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Tickle(LAXSubsequence):
    """
    Subsequence: Tickle

    Heat the ion by applying a tickle at the secular frequency.
    """
    name = 'tickle'

    def build_subsequence(self):
        self.setattr_device('tickle')

    def prepare_subsequence(self):
        self.time_tickle_mu = self.get_parameter('time_tickle_us', group='timing', override=False, conversion=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self):
        self.tickle.on()
        delay_mu(self.time_tickle_mu)
        self.tickle.off()
