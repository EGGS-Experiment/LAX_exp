import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.sequences import PMTDiscrimination


class pmt_disc_test(LAXExperiment, Experiment):
    """
    pmt_disc_test
    """

    name = "pmt_disc_test"

    def build_experiment(self):
        # timing
        self.setattr_argument("time_729_us",                        NumberValue(default=4000, ndecimals=5, step=1, min=1, max=10000000))

    def prepare_experiment(self):
        # get 729 beam
        self.setattr_device('qubit')


        # prepare sequences
        self.pmt_disc_sequence =                                    PMTDiscrimination(self)


    # MAIN SEQUENCE
    @kernel
    def run_main(self):
        self.core.reset()

        self.pmt_disc_sequence.run()
        self.core.break_realtime()

    def analyze(self):
        self.pmt_disc_sequence.analyze()
