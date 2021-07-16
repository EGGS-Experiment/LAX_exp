from artiq.experiment import *
import numpy as np
#import time

class Interferometer(EnvExperiment):
    """Record Interferometer Data"""

    def build(self):
        self.setattr_argument("num_samples", NumberValue(ndecimals=0, step=1))
        #make delay time and recordchannel arguments as well
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("suservo0")
        self.setattr_device("suservo0_ch0")

    def prepare(self):
        self.set_dataset("interferometer_data", np.full(num_samples, np.nan))

    @kernel
    def record(self):
        with self.core_dma.record("record"):
            for i in range(num_samples):
                self.mutate_dataset("interferometer_data", i, self.suservo0.get_adc(record_channel))
                delay(delay_time)

    @kernel
    def run(self):
        #initialize components
        self.core.reset()
        self.suservo0.init()

        #build record sequence
        self.record()
        record_handle = self.core_dma.get_handle("record")
        self.core.break_realtime()

        #record data
        self.core_dma.playback_handle(record_handle)