from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Linetrigger(LAXDevice):
    """
    Triggers an event off the AC line.
        Requires the AC line to already be converted into a TTL signal.
    """
    name = "linetrigger"

    # todo: have an ms to mu function instead maybe???
    parameters = {
        'gating_edge':                  ('linetrigger.gating_edge',                 None),
        'time_timeout_mu':              ('linetrigger.time_timeout_ms',             seconds_to_mu),
        'time_holdoff_mu':              ('linetrigger.time_holdoff_us',             seconds_to_mu),
    }
    core_devices = {
        'trigger': 'ttl4'
    }


    @kernel(flags={"fast-math"})
    def prepare_device(self):
        self.trigger.input()


    @kernel(flags={"fast-math"})
    def trigger(self):
        """
        Block until either the timeout period is reached, or a TTL edge is detected.
        """
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
