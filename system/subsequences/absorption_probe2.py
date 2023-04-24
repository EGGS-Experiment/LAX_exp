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
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')
        self.setattr_device('ttl0')

        self.setattr_argument("repetitions",                                NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000000))
        self.setattr_argument("repetitions_per_cooling",                    NumberValue(default=300, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_probe_us",                              NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_reset_us",                              NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

    def prepare_subsequence(self):
        self.time_doppler_cooling_mu =                                      self.get_parameter('time_doppler_cooling_us',
                                                                                               group='timing',
                                                                                               override=True,
                                                                                               conversion_function=seconds_to_mu, units=us)

        # tmp
        self._iter_loop_1 = np.array_split(np.arange(self.repetitions), self.repetitions / self.repetitions_per_cooling)
        self.counts_store = np.int32(0)

        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_reset_mu = self.core.seconds_to_mu(self.time_reset_us * us)

    @kernel(flags={"fast-math"})
    def run(self):

        # with parallel:
        # reset counts store
        self.counts_store = np.int32(0)

        # ensure 854nm on during sequence
        self.repump_qubit.on()

        # todo: move counter forward to schedule a very long pulse sequence
        # min_now = self.core.rtio_get_counter() + self.time_sequence_delay
        # if now_mu() < min_now:
        #     at_mu(min_now)

        # loop
        for _iter_arr in self._iter_loop_1:
            self.core.break_realtime()

            # DOPPLER COOLING
            # prepare beams for doppler cooling
            self.pump.cooling()
            self.repump_cooling.on()

            # doppler cool
            self.pump.on()
            delay_mu(self.time_doppler_cooling_mu)
            self.pump.off()

            # prepare probe beam
            self.pump.set_profile(1)

            # get counts
            for _iter_i in _iter_arr:

                # schedule pulse sequence
                self.repump_cooling.off()

                # apply probe beam
                self.pump.on()
                self.pmt.count(self.time_probe_mu)
                self.pump.off()

                # with parallel:
                # apply reset beam
                self.repump_cooling.on()
                delay_mu(self.time_reset_mu)

                self.counts_store += self.pmt.fetch_count()
