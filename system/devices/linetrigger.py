from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Linetrigger(LAXDevice):
    """
    Device: Linetrigger

    Triggers an event off the AC line.
    Requires the AC line to already be converted into a TTL signal.
    """
    name = "linetrigger"

    core_devices = {
        'trigger': 'ttl4'
    }

    def prepare_device(self):
        self.gating_edge = self.get_parameter('gating_edge', group='linetrigger', override=False)
        self.time_timeout_mu = self.get_parameter('time_timeout_ms', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=ms)
        self.time_holdoff_mu = self.get_parameter('time_holdoff_us', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # todo: tune delays
        self.trigger.input()

    @kernel(flags={"fast-math"})
    def trigger(self):
        """
        Block until either the timeout period is reached, or a TTL edge is detected.
        """
        # todo: tune delays
        # wait for PMT count
        self.trigger._set_sensitivity(1)
        time_input_mu = self.trigger.timestamp_mu(self.time_timeout_mu)

        # check if event has fired
        if time_input_mu > 0:

            # set RTIO time
            at_mu(time_input_mu)

            # close gating and add holdoff
            with parallel:
                self.trigger._set_sensitivity(0)
                delay_mu(self.holdoff_mu)

            return time_input_mu

        return -1
