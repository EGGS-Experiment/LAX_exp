import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.setattr_dataset("interferometer_data", np.full(10, np.nan))

    @kernel
    def run(self):
        #initialize devices
        self.core.reset()
        self.sampler0.init()

        print("h1")
        self.core.break_realtime()

        self.sampler0.set_gain_mu(0, 2)
        print("h2")
        self.core.break_realtime()

        print("h3")
        self.core.break_realtime()

        self.sampler0.sample(self.interferometer_data)

        print("h4")
        self.core.break_realtime()