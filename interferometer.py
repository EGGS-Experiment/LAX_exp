import numpy as np
from artiq.experiment import *

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        self.setattr_argument("num_samples", NumberValue(ndecimals=0, step=1))
        self.setattr_argument("delay_time", NumberValue(ndecimals=0, step=1))
        self.setattr_argument("record_channel", NumberValue(ndecimals=0, step=1))
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("sampler0")
        self.setattr_dataset("interferometer_data", np.full(100, 0))

    # @kernel
    # def record(self):
    #     with self.core_dma.record("record"):
    #         self.sampler0.sample(self.interferometer_data)
    #         for i in range(100):
    #             self.mutate_dataset("interferometer_data", i, self.sampler0.sample(self.record_channel))
    #             delay(self.delay_time * us)

    @kernel
    def run(self):
        #initialize devices
        self.core.reset()
        self.sampler0.init()
        self.sampler0.set_gain_mu(0, 2)
        self.core.break_realtime()

        self.sampler0.sample(self.interferometer_data)

        # #build record sequence
        # self.record()
        # record_handle = self.core_dma.get_handle("record")
        # self.core.break_realtime()
        #
        # #record data
        # self.core_dma.playback_handle(record_handle)