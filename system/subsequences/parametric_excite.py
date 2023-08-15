import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
# todo: argument for num counts


class ParametricExcite(LAXSubsequence):
    """
    Subsequence: Parametric Excite

    Apply parametric excitation to the ion and read out the timestamps.
    """
    name = 'parametric_excite'

    def build_subsequence(self):
        # number of input counts to listen for
        self.setattr_argument("num_counts",     NumberValue(default=10000, ndecimals=0, step=1, min=1, max=10000000), group=self.name)

        # get relevant devices
        self.setattr_device('pmt')
        self.setattr_device('trigger_rf')
        self.setattr_device('dds_modulation')


    def prepare_subsequence(self):
        # get input triggering parameters
        self.time_pmt_gating_mu =               self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu, units=us)
        self.time_rf_gating_mu =                self.get_parameter('time_rf_gating_ns', group='rf', override=False, conversion_function=seconds_to_mu, units=ns)
        self.time_rf_holdoff_mu =               self.get_parameter('time_rf_holdoff_us', group='rf', override=False, conversion_function=seconds_to_mu, units=us)

        # set additional holdoff for modulation
        self.time_mod_delay_mu =                self.core.seconds_to_mu(200 * us)

        # create holder array to store PMT counts
        self.timestamp_mu_list =                np.zeros(self.num_counts, dtype=np.int64)



    @kernel(flags={"fast-math"})
    def run(self) -> TArray(TInt64, 1):
        # trigger sequence off same phase of RF
        self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)

        # reset modulation DDS phase and turn modulation on
        self.dds_modulation.reset_phase()
        self.dds_modulation.on()
        # add extra delay for RTIO slack (min 200us), as well as to reduce effects of initial conditions
        delay_mu(self.time_mod_delay_mu)

        # get timestamped photon counts
        time_start_mu = now_mu()
        # timestamped_count_list = array(self.pmt.timestamp_counts(num_counts, self.time_pmt_gating_mu))
        self.pmt.timestamp_counts(self.timestamp_mu_list, self.time_pmt_gating_mu)

        # add slack and turn off modulation
        with parallel:
            self.core.break_realtime()
            self.dds_modulation.off()

        # remove initial time
        self.timestamp_mu_list -= time_start_mu

        # return timestamped counts
        return self.timestamp_mu_list
