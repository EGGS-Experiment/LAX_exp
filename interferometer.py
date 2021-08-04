import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.setattr_dataset("interferometer_data", np.full(10, 1))

    @kernel
    def run(self):
        #initialize devices
        self.core.reset()

        self.sampler0.init()
        self.sampler0.set_gain_mu(0, 2)

        print(self.interferometer_data)
        self.core.break_realtime()

        #sample
        self.sampler0.sample(self.interferometer_data)

        print("h4")
        self.core.break_realtime()