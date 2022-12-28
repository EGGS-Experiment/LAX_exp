from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class RFSync(LAXDevice):
    """
    Correlates PMT counts with the phase of the modulation signal.
    """
    name = "rf_sync"

    parameters = {
        'pmt_gating_edge':          ('pmt.pmt_gating_edge',                 None),
        'rf_gating_edge':           ('rf.rfsync_gating_edge',               None),

        'pmt_gating_timeout_mu':    ('pmt.pmt_gating_timeout_us',           seconds_to_mu),
        'rf_gating_timeout_mu':     ('rf.rfsync_gating_timeout_us',         seconds_to_mu),
    }
    core_devices = {
        'pmt':          'ttl0',
        'rf_sync':      'ttl3'
    }

    @kernel(flags={"fast-math"})
    def correlate_count(self, timeout_mu):
        """
        Block until either the timeout period is reached, or a TTL edge is detected.
        """
        # wait for PMT count
        time_end_pmt_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
        time_input_pmt_mu = self.pmt_counter.timestamp_mu(time_end_pmt_mu)

        # check if event has fired
        if time_input_pmt_mu > 0:

            # set RTIO time and add slack
            at_mu(time_input_pmt_mu)
            delay_mu(self.time_slack_mu)

            # get timestamp of RF event
            time_end_rf_mu = self.rf_sync.gate_rising_mu(self.time_timeout_rf_mu)
            time_input_rf_mu = self.rf_sync.timestamp_mu(time_end_rf_mu)

            # close input gating
            self.rf_sync.count(time_end_rf_mu)
            self.pmt_counter.count(time_end_pmt_mu)
            self.core.break_realtime()

            # add data to dataset
            self.core.break_realtime()
            return (time_end_pmt_mu, time_input_rf_mu)
