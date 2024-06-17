import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Ramsey, Readout, RescueIon


class RamseySpectroscopy(LAXExperiment, Experiment):
    """
    Experiment: Ramsey Spectroscopy

    Measures ion fluorescence after conducting a Ramsey Spectroscopy sequence.
    """
    name = 'Ramsey Spectroscopy'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=5, ndecimals=0, step=1, min=1, max=100000))

        # ramsey parameters
        self.setattr_argument("freq_ramsey_mhz_list",               Scannable(
                                                                        default=CenterScan(104.335, 0.5, 0.001),
                                                                        global_min=30, global_max=200, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ))
        # get devices
        self.setattr_device('qubit')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)
        self.ramsey_subsequence =                                   Ramsey(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # convert ramsey detunings to ftw
        self.freq_ramsey_ftw_list =                                 np.array([hz_to_ftw(freq_mhz * MHz)
                                                                              for freq_mhz in list(self.freq_ramsey_mhz_list)])

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_ramsey_ftw_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.ramsey_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # sweep ramsey detunings
            for freq_ftw in self.freq_ramsey_ftw_list:

                # set ramsey detuning
                self.qubit.set_mu(freq_ftw, asf=self.qubit.ampl_qubit_asf)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # do ramsey sequence
                self.ramsey_subsequence.run_dma()

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()
