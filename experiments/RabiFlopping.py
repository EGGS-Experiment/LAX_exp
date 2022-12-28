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
        # rabi flopping parameters
        self.setattr_argument("time_rabi_us_list",                  Scannable(
                                                                        default=LinearScan(0, 400, 10, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ))

        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=110, ndecimals=5, step=1, min=1, max=10000))

    def prepare_experiment(self):
        # rabi flopping timing
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [us_to_mu(max_time_us - time_us), us_to_mu(time_us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])

        # rabi flopping frequency
        self.freq_rabiflop_ftw =                                    mhz_to_ftw(self.freq_rabiflop_mhz)

        # get devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        # tmp remove
        self.rubbish_initialize =                               InitializeQubit(self)
        print(self.rubbish_initialize.dma_name)
        # tmp remove clear
        self.readout_subsequence =                                  Readout(self)

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(self.time_rabiflop_mu_list), 2)))
        self.setattr_dataset('results')


    # PREPARE MAIN SEQUENCE
    @kernel
    def run_initialize(self):
        self.core.reset()

        # set qubit beam parameters
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf)

        # record subsequences onto DMA
        # tmp remove
        self.rubbish_initialize.record_dma()
        # tmp remove clear
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        self.core.break_realtime()


    # MAIN LOOP
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
                self.qubit.cfg_sw(True)
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.cfg_sw(False)

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_pair_mu[1], self.pmt.fetch_count())
                    self.core.break_realtime()

    @rpc(flags={"async"})
    def update_dataset(self, time_mu, counts):
        self.results[self._result_iter] = np.array([time_mu, counts])
        self._result_iter += 1
        print(self._result_iter)
