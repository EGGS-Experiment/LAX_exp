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
                                                                        default=RangeScan(0, 400, 1001, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ))

        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=110, ndecimals=5, step=1, min=1, max=10000))

    def prepare_experiment(self):
        # rabi flopping timing
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [seconds_to_mu((max_time_us - time_us) * us), seconds_to_mu(time_us * us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])

        # rabi flopping frequency
        self.freq_rabiflop_mhz =                                    mhz_to_ftw(self.freq_rabiflop_mhz)

        # get devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(self.time_rabiflop_mu_list), 2)))
        self.setattr_dataset('results')


    # PREPARE MAIN SEQUENCE
    @kernel
    def run_initialize(self):
        self.core.reset()

        # set qubit beam parameters
        self.qubit.set_mu(self.freq_rabiflop_mhz, asf=self.qubit.ampl_qubit_asf)

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.core.break_realtime()

        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # load subsequences from DMA
        self.initialize_subsequence.load_dma()
        self.core.break_realtime()

        self.readout_subsequence.load_dma()
        self.core.break_realtime()


    # MAIN LOOP
    @kernel
    def loop(self):
        # sweep time
        for time_rabi_pair_mu in self.time_rabiflop_mu_list:

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
                self.update_dataset(time_rabi_pair_mu[1], self.pmt_counter.fetch_count())
                self.core.break_realtime()
