import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")

    @kernel
    def run(self):
        #initialize
        self.core.reset()
        self.set_dataset("interferometer_data", np.full(10, np.nan), broadcast=True)

        self.sampler0.init()
        self.sampler0.set_gain_mu(0, 2)

        self.core.break_realtime()
        holder = [0]*8

        for ind in range(100):
            self.sampler0.sample(holder)
            self.core.mutate_dataset("interferometer_data", ind, holder[0])
            delay(100*us)