from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class TriggerLine(LAXDevice):
    """
    Device: Trigger Line

    Triggers an event off the AC line.
    Requires the AC line to already be converted into a TTL signal.
    """
    name = "trigger_line"
    core_device = ('ttl_input', 'ttl5')
    kernel_invariants = {
        "time_timeout_mu",
        "time_holdoff_mu"
    }

    def prepare_device(self):
        # get triggering parameters
        self.time_timeout_mu = self.get_parameter('time_timeout_ms', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=ms)
        self.time_holdoff_mu = self.get_parameter('time_holdoff_us', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def trigger(self, time_gating_mu: TInt64, time_holdoff_mu: TInt64) -> TInt64:
        """
        Trigger off a rising edge of the AC line.
        Times out if no edges are detected.
        Arguments:
            time_gating_mu  (int)   : the maximum waiting time (in machine units) for the trigger signal.
            time_holdoff_mu (int64) : the holdoff time (in machine units)
        Returns:
                            (int64) : the input time of the trigger signal.
        """
        # wait for line trigger input
        time_trigger_mu = self.ttl_input.timestamp_mu(self.ttl_input.gate_rising_mu(time_gating_mu))

        # ensure input timestamp is valid
        if time_trigger_mu > 0:
            # set rtio hardware time to input trigger time
            at_mu(time_trigger_mu + time_holdoff_mu)
            return time_trigger_mu

        # reset RTIO if we don't get receive trigger signal for some reason
        else:
            # add slack before resetting system
            self.core.break_realtime()
            self.ttl_input._set_sensitivity(0)
            self.core.reset()

        # return -1 if we time out
        return -1
