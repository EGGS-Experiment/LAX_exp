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
        # self.setattr_device('pump')
        self.setattr_device('trigger_rf')
        self.setattr_device('dds_modulation')


    def prepare_subsequence(self):
        # get input triggering parameters
        self.time_pmt_gating_mu =               self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu, units=us)
        self.time_rf_gating_mu =                self.get_parameter('time_rf_gating_ns', group='rf', override=False, conversion_function=seconds_to_mu, units=ns)
        self.time_rf_holdoff_mu =               self.get_parameter('time_rf_holdoff_us', group='rf', override=False, conversion_function=seconds_to_mu, units=us)

        # set additional holdoff for modulation
        self.time_mod_delay_mu =                self.core.seconds_to_mu(200 * us)


    @kernel(flags={"fast-math"})
    def run(self, num_counts: TInt32) -> TArray(TInt64, 1):
        # trigger sequence off same phase of RF
        self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)

        # reset modulation DDS phase and turn modulation on
        self.dds_modulation.reset_phase()
        self.dds_modulation.on()

        # add extra delay for slack (min 200us) and to reduce effects of initial conditions
        delay_mu(self.time_mod_delay_mu)

        # get timestamped photon counts
        time_start_mu = now_mu()
        timestamped_count_list = array(self.pmt.timestamp_counts(num_counts, self.time_pmt_gating_mu))

        # add slack and turn off modulation
        self.core.break_realtime()
        self.dds_modulation.off()

        # remove initial time
        timestamped_count_list -= time_start_mu

        # return timestamped counts
        return timestamped_count_list

    # @kernel(flags={"fast-math"})
    # def run(self, num_counts: TInt32) -> TArray(TInt64, 1):
    #     # trigger sequence off same phase of RF
    #     time_tmp_mu = self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)
    #
    #     self.dds_modulation.on()
    #
    #     # reset modulation DDS phase and turn modulation on
    #     self.dds_modulation.reset_phase()
    #
    #     # add extra delay for slack (min 200us) and to reduce effects of initial conditions
    #     delay_mu(self.time_mod_delay_mu)
    #
    #     # get timestamped photon counts
    #     time_start_mu = now_mu()
    #     timestamped_count_list = array(self.pmt.timestamp_counts(num_counts, self.time_pmt_gating_mu))
    #
    #     # add slack and turn off modulation
    #     self.core.break_realtime()
    #     self.dds_modulation.off()
    #
    #     # remove initial time
    #     timestamped_count_list -= time_start_mu
    #
    #     # return timestamped counts
    #     return timestamped_count_list
