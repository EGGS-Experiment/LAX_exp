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
        self.set_dataset("interferometer_data", np.full(10, np.nan), broadcast=False)

        self.core.break_realtime()

        holder = [0]*8
        self.sampler0.init()

        for i in range(8):
            self.sampler0.set_gain_mu(7-i, 0)

        self.core.break_realtime()

        for ind in range(10):
            self.sampler0.sample_mu(holder)
            self.mutate_dataset("interferometer_data", ind, holder[0])
            delay(56*us)

        print("done")