from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class PMTProcessData(EnvExperiment):
    """
    Print dataset of PMT.
    """

    def build(self):
        self.setattr_device("core")

    def prepare(self):
        #self.time_step_mu = self.core.seconds_to_mu(1 * us)
        # set dataset
        self.dataset_name = 'oktc'
        #self.set_dataset(self.dataset_name, np.zeros(self.num_bins), broadcast=True

    @kernel
    def run(self):
        pass

    def analyze(self):
        th1 = self.get_dataset(self.dataset_name)
        print(th1)
        #np.savetxt('oktc.csv', th1, delimiter=',')
