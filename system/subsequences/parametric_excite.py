from numpy import array
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
        # get relevant devices
        self.setattr_device('pmt')
        self.setattr_device('pump')
        self.setattr_device('trigger_rf')
        self.setattr_device('dds_modulation')


    def prepare_subsequence(self):
        # get input triggering parameters
        self.time_pmt_gating_mu =               self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu, units=us)
        self.time_rf_gating_mu =                self.get_parameter('time_rf_gating_ns', group='rf', override=False, conversion_function=seconds_to_mu, units=ns)
        self.time_rf_holdoff_mu =               self.get_parameter('time_rf_holdoff_us', group='rf', override=False, conversion_function=seconds_to_mu, units=us)

        # set additional holdoff for modulation
        self.time_mod_delay_mu =                self.core.seconds_to_mu(100 * us)


    @kernel(flags={"fast-math"})
    def run(self, num_counts: TInt32) -> TArray(TInt64, 1):
        # turn on cooling beam and allow ion to recool
        self.pump.on()
        delay_mu(self.time_cooling_holdoff_mu)

        # trigger sequence off same phase of RF
        self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)
        at_mu(now_mu())

        # reset modulation DDS phase
        self.mod_dds.cpld.io_update.pulse_mu(8)
        # turn modulation on, then wait a given time to reduce effect of initial conditions
        self.mod_dds.on()
        delay_mu(self.time_mod_delay_mu)

        # get timestamped photon counts
        time_start_mu = now_mu()
        timestamped_count_list = self.pmt.timestamp_counts(num_counts, self.time_pmt_gating_mu)

        # turn off modulation and reset
        self.mod_dds.cfg_sw(False)
        self.core.reset()

        # return timestamped counts
        return array(timestamped_count_list) - time_start_mu
