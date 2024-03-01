from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class TriggerRF(LAXDevice):
    """
    Device: Trigger RF

    Wrapper for the ttl_input object that reads in the RF synchronization signal.
    """
    name = "trigger_rf"
    core_device = ('ttl_input', 'ttl7')
    kernel_invariants = set()

    def prepare_device(self):
        pass

    @kernel(flags={"fast-math"})
    def trigger(self, time_gating_mu: TInt64, time_holdoff_mu: TInt64) -> TInt64:
        """
        Trigger off a rising edge of the RF sync input.
        Times out if no edges are detected.
        Arguments:
            time_gating_mu  (int)   : the maximum waiting time (in machine units) for the trigger signal.
            time_holdoff_mu (int64) : the holdoff time (in machine units)
        Returns:
                            (int64) : the input time of the trigger signal.
        """
        # trigger sequence off same phase of RF
        self.ttl_input.gate_rising_mu(time_gating_mu)
        time_trigger_mu = self.ttl_input.timestamp_mu(now_mu())

        # ensure input timestamp is valid
        if time_trigger_mu >= 0:
            # activate modulation and enable photon counting
            at_mu(time_trigger_mu + time_holdoff_mu)
            return now_mu()

        # reset RTIO if we don't get receive trigger signal for some reason
        else:
            # add slack before resetting system
            self.core.break_realtime()
            self.ttl_input._set_sensitivity(0)
            self.core.reset()

        # return -1 if we time out
        return -1
