import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """
    Interferometer Run
    Records the data from the interferometer.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.set_dataset("interferometer_data", np.zeros(100), broadcast=True)

    @kernel
    def run(self):
        #initialize
        self.core.reset()


        self.core.reset()

        holder = [0]*8
        self.sampler0.init()
        self.core.reset()
        for i in range(8):
            self.sampler0.set_gain_mu(7-i, 0)
            self.core.break_realtime()

        self.core.reset()

        for ind in range(100):
            self.sampler0.sample_mu(holder)
            self.mutate_dataset("interferometer_data", ind, holder[0])
            self.core.reset()

        #print("done")