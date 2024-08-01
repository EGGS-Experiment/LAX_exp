import extensions
import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
import labrad

class PMTTest(LAXExperiment, Experiment):
    """
    Experiment: PMT Test

    Test PMT
    """
    name = 'PMT Test'
    DATASET_KEY = 'test_pmt_counts'

    def build_experiment(self):

        self.setattr_device('pmt')

    def prepare_experiment(self):
        self.set_dataset(self.DATASET_KEY, [])

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        pass


    @kernel(flags={"fast-math"})
    def run_main(self):

        for i in range(100):
            self.pmt.count(self.pmt_sample_time_us)
            self.append_to_dataset(self.DATASET_KEY, self.pmt.fetch_count())
            self.core.break_realtime()

    # ANALYSIS
    def analyze_experiment(self):

        counts = self.get_dataset(self.DATASET_KEY)
        print(counts)


