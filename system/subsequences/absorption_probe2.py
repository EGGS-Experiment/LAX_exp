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
    name = 'absorption_probe'

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

        self.setattr_argument("repetitions_per_point",                      NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000000))
        self.setattr_argument("repetitions_per_cooling",                    NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_probe_us",                              NumberValue(default=2, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_reset_us",                              NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

    def prepare_subsequence(self):
        self.time_doppler_cooling_mu =      self.get_parameter('time_doppler_cooling_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_probe_mu =                self.get_parameter('time_probe_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # tmp
        self._iter_loop_1 = np.arange(self.repetitions_per_point)
        self._iter_loop_2 = np.arange(self.repetitions_per_cooling)

    @kernel(flags={"fast-math"})
    def run(self):

        counts_store = np.int32(0)

        # loop
        for i in self._iter_loop_1:

            # DOPPLER COOLING
            # prepare beams for doppler cooling
            self.repump_cooling.on()
            self.pump.cooling()

            # doppler cool
            self.pump.on()
            delay_mu(self.time_doppler_cooling_mu)
            self.pump.off()

            # prepare probe beam
            self.pump.set_profile(1)
            self.repump_cooling.off()


            # GET PROBE FLUORESCENCE
            # get probe fluorescence
            for j in self._iter_loop_2:

                # record probe fluorescence
                self.pump.on()
                self.pmt.count(self.time_probe_mu)
                self.pump.off()

                # repump ion
                with parallel:
                    with sequential:
                        self.repump_cooling.on()
                        delay_mu(self.time_reset_us)
                        self.repump_cooling.off()
                    counts_store += self.pmt.fetch_count()

        return float(counts_store) / float(self.repetitions_per_point)
