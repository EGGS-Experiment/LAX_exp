import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon


class RabiFlopping(LAXExperiment, Experiment):
    """
    Experiment: Rabi Flopping

    Measures ion fluorescence vs 729nm pulse time and frequency.
    """

    name = 'Rabi Flopping'

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=30, ndecimals=0, step=1, min=1, max=10000))

        # rabi flopping arguments
        self.setattr_argument("time_rabi_us_list",                  Scannable(
                                                                        default=RangeScan(0, 250, 251, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ), group=self.name)
        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=104.335, ndecimals=5, step=1, min=1, max=10000), group=self.name)

        # get devices
        self.setattr_device('qubit')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # convert rabi flopping arguments
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [seconds_to_mu((max_time_us - time_us) * us), seconds_to_mu(time_us * us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])
        self.freq_rabiflop_ftw =                                    hz_to_ftw(self.freq_rabiflop_mhz * MHz)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.reset()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # set qubit beam parameters
        #self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf)
        # todo: check if this gives us problems like before related to profile=0
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
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
                    self.update_results(time_rabi_pair_mu[1], self.readout_subsequence.fetch_count())
                    self.core.break_realtime()
