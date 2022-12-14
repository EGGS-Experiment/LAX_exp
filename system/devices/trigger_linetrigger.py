from artiq.experiment import *

from LAX_exp.utilities.conversions import *
from LAX_exp.base import LAXDevice


class Linetrigger(LAXDevice):
    """
    Triggers an event off the AC line.
        Requires the AC line to already be converted into a TTL signal.
    """
    name = "linetrigger"

    parameters = {
        'gating_edge':              ('triggers.linetrigger_gating_edge',            None),
        'gating_timeout_mu':        ('triggers.linetrigger_gating_timeout_mu',      seconds_to_mu),
    }
    core_devices = {
        'trigger': 'ttl4'
    }


    @kernel(flags='fast-math')
    def prepare_device(self):
        self.core.break_realtime()
        self.trigger.input()


    @kernel(flags='fast-math')
    def trigger(self):
        """
        Block until either the timeout period is reached, or a TTL edge is detected.
        """
        # wait for PMT count
        time_end_mu = self.pmt_counter.gate_rising_mu(self.gating_timeout_mu)
        time_input_mu = self.pmt_counter.timestamp_mu(time_end_mu)

        # check if event has fired
        if time_input_mu > 0:

            # set RTIO time and add slack
            at_mu(time_input_mu)
            delay_mu(5000)

            # wait until end of period
            self.rf_sync.count(time_end_mu)
            return (True, time_input_mu)

        return (False, time_input_mu)

    # todo: check this decorator is correct
    @rpc
    def set_trigger_timeout_ms(self, timeout_ms):
        """
        Change the trigger timeout.
        """
        self.gating_timeout_mu = seconds_to_mu(timeout_ms * ms)
