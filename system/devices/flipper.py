from artiq.experiment import *

from LAX_exp.base import LAXDevice
from LAX_exp.extensions import *

from os import environ
import labrad
from artiq.language.units import *


class Flipper(LAXDevice):
    """
    High-level api functions for using the flipper
    """

    name = 'flipper'
    core_device = ('flipper', 'ttl15')

    def prepare_device(self):

        self.wait_time_mu = us_to_mu(1000e3)
        self.time_flipper_trigger_mu = self.core.seconds_to_mu(50*ms)


    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @kernel
    def flip(self):
        self.core.break_realtime()
        self.flipper.on()
        delay_mu(self.time_flipper_trigger_mu)
        self.flipper.off()
        self.core.wait_until_mu(now_mu())
        delay_mu(self.wait_time_mu)
