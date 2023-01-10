import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.sequences import PMTDiscrimination


class pmt_disc_test(LAXExperiment, Experiment):
    """
    pmt_disc_test
    """

    name = 'pmt_disc_test"

    def build_experiment(self):
        # timing
        self.setattr_argument("time_729_us",                        NumberValue(default=4000, ndecimals=5, step=1, min=1, max=10000000))

    def prepare_experiment(self):
        # get 729 beam
        self.setattr_device('qubit')

        # convert frequencies to machine units
        self.freq_qubit_scan_ftw =                                  np.array([mhz_to_ftw(freq_mhz) for freq_mhz in self.freq_qubit_scan_mhz])

        # prepare subsequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.rabiflop_subsequence =                                 RabiFlop(self, time_rabiflop_us=self.time_729_us)
        self.readout_subsequence =                                  Readout(self)

        # prepare sequences
        self.pmt_disc_sequence =                                    PMTDiscrimination(self)

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(list(self.freq_qubit_scan_mhz)), 2)))
        self.setattr_dataset('results')


    # MAIN SEQUENCE
    @kernel
    def run_main(self):
        self.core.reset()

        self.pmt_disc_sequence.run()
        self.core.break_realtime()

        self.
