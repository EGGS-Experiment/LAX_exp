import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class ReadoutBackgroundSubtract(LAXSubsequence):
    """
    Subsequence: Readout Background Subtract

    Read out the background-subtracted ion fluorescence by getting counts with both
    the pump and repump beams on, then getting couts with only the pump on.
    """
    name = 'readout_background_subtract'
    kernel_invariants = {
        "time_readout_mu"
    }

    def build_subsequence(self):
        # sampling
        self.setattr_argument('time_sample_us',                 NumberValue(default=3000, ndecimals=1, step=500, min=100, max=100000), group='sampling')
        self.setattr_argument('signal_samples_per_point',       NumberValue(default=20, ndecimals=0, step=10, min=1, max=100), group='sampling')
        self.setattr_argument('background_samples_per_point',   NumberValue(default=5, ndecimals=0, step=2, min=1, max=100), group='sampling')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        # convert relevant timings to machine units
        self.time_slack_us =        2
        self.time_per_point_s =     (self.time_slack_us * us + self.time_sample_us * us)* (self.signal_samples_per_point + self.background_samples_per_point)

        self.time_slack_mu =        self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_sample_mu =       self.core.seconds_to_mu(self.time_sample_us * us)
        self.time_per_point_mu =    self.core.seconds_to_mu(self.time_per_point_s)

        # calculate the number of repetitions
        self.repetitions =          round(self.time_total_s / self.core.mu_to_seconds(self.time_per_point_mu))

        # create holder variable that can store averaged counts
        self._counts_signal =       np.int32(0)
        self._counts_background =   np.int32(0)

        # declare the loop iterators ahead of time to reduce overhead
        self._iter_repetitions =    np.arange(self.repetitions)
        self._iter_signal =         np.arange(self.signal_samples_per_point)
        self._iter_background =     np.arange(self.background_samples_per_point)

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform
        self.pump.readout()

        # readout pulse
        self.pump.on()
        self.repump_cooling.on()
        self.pmt.count(self.time_readout_mu)
        self.pump.off()

    @kernel(flags={"fast-math"})
    def readout_signal(self):
        # record alignment sequence - signal
        with self.core_dma.record('_PMT_ALIGNMENT_SIGNAL'):
            # activate doppler repump beam
            self.repump_cooling.on()

            # get signal counts
            for i in self._iter_signal:
                self.pmt.count(self.time_sample_mu)
                delay_mu(self.time_slack_mu)

        # record alignment sequence - background
        with self.core_dma.record('_PMT_ALIGNMENT_BACKGROUND'):
            # disable doppler repump beam
            self.repump_cooling.off()

            # get background counts
            for i in self._iter_background:
                self.pmt.count(self.time_sample_mu)
                delay_mu(self.time_slack_mu)

    @kernel(flags={"fast-math"})
    def fetch_counts_signal(self) -> TInt32:
        """
        Convenience function so that users don't have to separately instantiate the PMT
        device object to read counts.
        """
        return self.pmt.fetch_count()

    @kernel(flags={"fast-math"})
    def fetch_counts_background(self) -> TInt32:
        """
        Convenience function so that users don't have to separately instantiate the PMT
        device object to read counts.
        """
        return self.pmt.fetch_count()
