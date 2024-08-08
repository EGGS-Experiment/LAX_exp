from artiq.experiment import *

from LAX_exp.base import LAXDevice

from os import environ
import labrad


class Flipper(LAXDevice):
    """
    High-level api functions for using the flipper
    """

    name = "flipper"

    def prepare_device(self):
        core_device = ('flipper', 'ttl15')

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @kernel
    def flip(self):
        self.flipper.off()
        delay_mu(1000000)
        self.flipper.pulse_mu(10000000)
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
