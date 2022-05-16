import numpy as np
from artiq.experiment import *
import time

class Testing(EnvExperiment):
    """Testing"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl4")

    def prepare(self):
        self.set_dataset("temptest", np.full(10, np.nan), broadcast=False)
        self.dataset_name = "temptest"
        self.time_measure_mu = self.core.seconds_to_mu(50 * us)

    @kernel
    def record(self):
        self.core.break_realtime()
        with self.core_dma.record("tmp_exp"):
            self.ttl0.gate_rising_mu(self.time_measure_mu)
            #self.ttl4.off()

    @kernel
    def run(self):

    @kernel
    def run(self):
        self.core.reset()
        # record sequence
        self.record()
        self.core.break_realtime()
        # get handle
        handle = self.core_dma.get_handle('tmp_exp')
        for i in range(10):
            self.mutate_dataset("interferometer_data", i, self.ttl0.count)
            delay(100 * us)
