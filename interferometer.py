import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        # self.setattr_argument("num_samples", NumberValue(ndecimals=0, step=1))
        # self.setattr_argument("delay_time", NumberValue(ndecimals=0, step=1))
        # self.setattr_argument("record_channel", NumberValue(ndecimals=0, step=1))
        self.setattr_device("sampler0")
        self.setattr_dataset("interferometer_data", np.full(10, np.nan))

    @kernel
    def run(self):
        #initialize devices
        self.core.reset()
        self.sampler0.init()
        self.sampler0.set_gain_mu(0, 2)
        self.core.break_realtime()

        self.sampler0.sample(self.interferometer_data)
