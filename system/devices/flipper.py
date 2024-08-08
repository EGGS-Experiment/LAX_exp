from artiq.experiment import *

from LAX_exp.base import LAXDevice
from LAX_exp.extensions import *

from os import environ
import labrad


class Flipper(LAXDevice):
    """
    High-level api functions for using the flipper
    """

    name = 'flipper'
    core_device = ('flipper', 'ttl15')

    def prepare_device(self):

        self.wait_time_mu = us_to_mu(200e3)


    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @kernel
    def flip(self):
        self.core.break_realtime()
        self.flipper.off()
        delay_mu(1000000)
        self.flipper.pulse_mu(10000000)
        delay_mu(self.wait_time_mu)
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

