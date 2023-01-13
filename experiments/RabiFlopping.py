import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout


class RabiFlopping2(LAXExperiment, Experiment):
    """
    Rabi Flopping 2
    Measures ion fluorescence vs 729nm pulse time and frequency.
    """

    name = 'Rabi Flopping 2'

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # rabi flopping parameters
        self.setattr_argument("time_rabi_us_list",                  Scannable(
                                                                        default=RangeScan(0, 400, 401, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ))

        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=104.335, ndecimals=5, step=1, min=1, max=10000))

    def prepare_experiment(self):
        # get devices
        self.setattr_device('qubit')

        # rabi flopping timing
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [us_to_mu(max_time_us - time_us), us_to_mu(time_us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])

        # rabi flopping frequency
        self.freq_rabiflop_ftw =                                    mhz_to_ftw(self.freq_rabiflop_mhz)

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(self.time_rabiflop_mu_list), 2)))
        self.setattr_dataset('results')


    # MAIN SEQUENCE
    @kernel
    def run_initialize(self):
        self.core.reset()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # set qubit beam parameters
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)
        self.core.break_realtime()

    @kernel
    def run_main(self):
        for trial_num in range(self.repetitions):

            # sweep time
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # wait given time
                delay_mu(time_rabi_pair_mu[0])

                # rabi flopping w/qubit laser
                self.qubit.on()
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.off()

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_pair_mu[1], self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

    @rpc(flags={"async"})
    def update_dataset(self, time_mu, counts):
        self.results[self._result_iter] = np.array([self.core.mu_to_seconds(time_mu), counts])
        self._result_iter += 1
