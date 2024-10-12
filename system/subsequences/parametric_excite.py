import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class ParametricExcite(LAXSubsequence):
    """
    Subsequence: Parametric Excite

    Apply parametric excitation to the ion and read out the timestamps.
    """
    name = 'parametric_excite'
    kernel_invariants = {
        "time_pmt_gating_mu",
        "time_rf_gating_mu",
        "time_rf_holdoff_mu"
    }

    def build_subsequence(self):
        # number of input counts to listen for
        self.setattr_argument("num_counts",     NumberValue(default=5000, precision=0, step=1, min=1, max=10000000), group=self.name)

        # get relevant devices
        self.setattr_device('pmt')
        self.setattr_device('trigger_rf')
        self.setattr_device('dds_parametric')

    def prepare_subsequence(self):
        # get input triggering parameters
        self.time_pmt_gating_mu =   self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu, units=us)
        self.time_rf_gating_mu =    self.get_parameter('time_rf_gating_ns', group='rf', override=False, conversion_function=seconds_to_mu, units=ns)
        self.time_rf_holdoff_mu =   self.get_parameter('time_rf_holdoff_us', group='rf', override=False, conversion_function=seconds_to_mu, units=us)

        # create holder array to store PMT counts
        self.timestamp_mu_list =    np.zeros(self.num_counts, dtype=np.int64)

    @kernel(flags={"fast-math"})
    def run(self) -> TArray(TInt64, 1):
        # trigger sequence off same phase of RF
        self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)

        # reset modulation DDS phase and turn modulation on
        self.dds_parametric.reset_phase()
        self.dds_parametric.on()

        # get timestamped photon counts
        self.pmt.timestamp_counts(self.timestamp_mu_list, self.time_pmt_gating_mu)

        # add slack and turn off modulation
        self.core.break_realtime()
        self.dds_parametric.off()

        # return timestamped counts
        return self.timestamp_mu_list
