from artiq.experiment import *
from os import environ
import labrad
from artiq.language.units import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


import labrad
from os import environ

class Flipper(LAXDevice):
    """
    High-level API functions for using the flipper
    """
    name = 'flipper'
    core_device = ('flipper', 'ttl15')

    def prepare_device(self):
        # define times so flipper finishes before other code is called
        self.wait_time_mu = self.core.seconds_to_mu(1000 * ms)
        self.time_flipper_trigger_mu = self.core.seconds_to_mu(10 * ms)

    @kernel(flags={"fast-math"})
    def flip(self) -> TNone:
        """
        Flip the mirror so light is sent to either PMT or camera
        """
        # ensure flipper TTL is initialized to 0
        self.flipper.off()
        delay_mu(self.time_flipper_trigger_mu)

        # pulse flipper TTL
        self.flipper.pulse_mu(self.time_flipper_trigger_mu)

        # add delay to allow flipper device to latch
        delay_mu(self.wait_time_mu)
        # define times so flipper finishes before other code is called
        self.wait_time_mu = self.core.seconds_to_mu(1000 * ms)
        self.time_flipper_trigger_mu = self.core.seconds_to_mu(10 * ms)

    @kernel(flags={"fast-math"})
    def flip(self) -> TNone:
        """
        Flip the mirror so light is sent to either PMT or camera
        """
        # ensure flipper TTL is initialized to 0
        self.flipper.off()
        delay_mu(self.time_flipper_trigger_mu)

        # pulse flipper TTL
        self.flipper.pulse_mu(self.time_flipper_trigger_mu)

        # add delay to allow flipper device to latch
        delay_mu(self.wait_time_mu)

        # synchronize with timeline
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()