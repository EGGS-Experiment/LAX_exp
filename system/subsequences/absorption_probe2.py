from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

# tmp remove
import numpy as np
# tmp remove


class AbsorptionProbe2(LAXSubsequence):
    """
    Subsequence: Absorption Probe2

    Cool the ion before switching to a short, low-power pulse to measure absorption on a transition.
    """
    name = 'absorption_probe2'
    kernel_invariants = {
        "time_doppler_cooling_mu",
        "time_probe_mu",
        "time_reset_mu",
        "_loop_iter",
    }

    def build_subsequence(self):
        # subsequence arguments
        self.setattr_argument("repetitions_per_point",  NumberValue(default=50, precision=0, step=1, min=1, max=1000), group='absorption_probe')
        self.setattr_argument("time_probe_us",          NumberValue(default=2, precision=0, step=1, min=1, max=10000), group='absorption_probe')
        self.setattr_argument("time_reset_us",          NumberValue(default=10, precision=0, step=1, min=1, max=10000), group='absorption_probe')

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        # get doppler cooling time
        self.time_doppler_cooling_mu =  self.get_parameter('time_doppler_cooling_us',
                                                           group='timing',
                                                           override=False,
                                                           conversion_function=seconds_to_mu, units=us)

        # convert pulse times to mu
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_reset_mu = self.core.seconds_to_mu(self.time_reset_us * us)

        # initialize loop variables
        self.counts_store = np.int32(0)
        self._loop_iter =   np.arange(self.repetitions_per_point)


    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare probe beam
        self.pump.set_profile(1)

        # get counts for multiple loops
        for i in self._loop_iter:

            # turn repump off
            self.repump_cooling.off()

            # apply probe beam and read counts
            self.pump.on()
            self.pmt.count(self.time_probe_mu)
            self.pump.off()

            # apply reset beam
            self.repump_cooling.on()
            delay_mu(self.time_reset_mu)
            self.repump_cooling.off()

    @kernel(flags={"fast-math"})
    def get_counts(self) -> TInt32:
        """
        Retrieve stored counts from memory.
        """
        # reset counts store
        self.counts_store = np.int32(0)
        self.core.break_realtime()

        # retrieve PMT counts and combine
        for i in self._loop_iter:
            self.counts_store += self.pmt.fetch_count()

        return self.counts_store
