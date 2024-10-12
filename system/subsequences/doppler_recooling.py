import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class DopplerRecooling(LAXSubsequence):
    """
    Subsequence: Doppler Recooling

    Apply parametric excitation to the ion and read out the timestamps.
    """
    name = 'doppler_recooling'
    kernel_invariants = {
        "time_readout_mu",
        "time_pmt_gating_mu"
    }

    def build_subsequence(self):
        # number of input counts to listen for
        self.setattr_argument("num_counts",     NumberValue(default=200, precision=0, step=1, min=10, max=10000000), group=self.name)

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        # get timing parameters
        self.time_readout_mu =      self.get_parameter('time_readout_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_pmt_gating_mu =   self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu, units=us)

        # create holder array to store PMT counts
        self.timestamp_mu_list =    np.zeros(self.num_counts, dtype=np.int64)

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform
        self.pump.readout()
        self.pump.on()
        self.repump_cooling.on()

        # get start time
        time_start_mu = now_mu()

        # get timestamped photon counts
        self.pmt.timestamp_counts(self.timestamp_mu_list, self.time_pmt_gating_mu)

        # ensure readout beam is on for at least the minimum time
        at_mu(time_start_mu + self.time_readout_mu)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def fetch_count(self) -> TArray(TInt64, 1):
        """
        Convenience function to decouple the run method from result retrieval.
        """
        return self.timestamp_mu_list
