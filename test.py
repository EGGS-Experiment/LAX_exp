import numpy as np
from artiq.experiment import *
import time

class Interferometer(EnvExperiment):
    """Testing"""

    def build(self):
        # self.setattr_argument("num_samples", NumberValue(ndecimals=0, step=1))
        # self.setattr_argument("delay_time", NumberValue(ndecimals=0, step=1))
        # self.setattr_argument("record_channel", NumberValue(ndecimals=0, step=1))
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("suservo0")
        self.setattr_device("suservo0_ch0")
        #self.setattr_dataset("interferometer_data", np.full(self.num_samples, np.nan))
        self.set_dataset("interferometer_data", np.full(100, np.nan), broadcast=True)

    @kernel
    def record(self):
        with self.core_dma.record("record"):
            for i in range(100):
                self.mutate_dataset("interferometer_data", i, self.suservo0.get_adc(0))
                delay(51 * us)

    @kernel
    def run(self):
        #initialize devices
        self.core.reset()
        self.suservo0.init()
        # self.suservo0.set_config(1)
        self.suservo0.set_pgia_mu(0, 0)
        self.core.break_realtime()

        #build record sequence
        self.record()
        record_handle = self.core_dma.get_handle("record")
        self.core.break_realtime()

        #record data
        self.core_dma.playback_handle(record_handle)
