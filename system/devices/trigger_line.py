from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class TriggerLine(LAXDevice):
    """
    Device: Line Trigger

    Triggers an event off the AC line.
    Requires the AC line to already be converted into a TTL signal.
    """
    name = "trigger_line"
    core_device = ('trigger', 'ttl5')

    def prepare_device(self):
        # get triggering parameters
        self.time_timeout_mu = self.get_parameter('time_timeout_ms', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=ms)
        self.time_holdoff_mu = self.get_parameter('time_holdoff_us', group='linetrigger', override=False, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def trigger(self):
    # def trigger(self, time_gating_mu: TInt64, time_holdoff_mu: TInt64) -> TInt64:
        """
        Block until either the timeout period is reached, or a TTL edge is detected.
        """
        self.core.break_realtime()

        # wait for trigger input
        time_trigger_mu = self.counter_trigger.timestamp_mu(self.trigger.gate_rising_mu(self.time_timeout_mu))

        # respond to trigger input
        if time_trigger_mu > 0:

            # set rtio hardware time to input trigger time
            at_mu(time_trigger_mu + self.time_holdoff_mu)
            return time_trigger_mu

        # if no input detected, return -1
        return -1
        # # trigger sequence off same phase of RF
        # self.input.gate_rising_mu(time_gating_mu)
        # time_trigger_mu = self.input.timestamp_mu(now_mu())
        #
        # # ensure input timestamp is valid
        # if time_trigger_mu >= 0:
        #
        #     # activate modulation and enable photon counting
        #     at_mu(time_trigger_mu + time_holdoff_mu)
        #     return now_mu()
        #
        # # reset RTIO if we don't get receive trigger signal for some reason
        # else:
        #     # add holdoff delay before resetting system
        #     self.core.break_realtime()
        #     self.rf_clock._set_sensitivity(0)
        #     self.core.reset()
        #
        # return -1
