from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.sequences import PMTDiscrimination


class pmt_disc_test(LAXExperiment, Experiment):
    """
    pmt_disc_test
    """

    name = "pmt_disc_test"

    def prepare_experiment(self):
        self.pmt_disc_sequence = PMTDiscrimination(self)


    # MAIN SEQUENCE
    @kernel
    def run_main(self):
        self.core.reset()
        self.pmt_disc_sequence.run()
        self.core.break_realtime()

    def analyze(self):
        self.pmt_disc_sequence.analyze()
