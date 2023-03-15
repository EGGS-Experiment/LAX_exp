from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class DACFast(EnvExperiment):
    """
    DACFast
    Set a value on the Fastino.
    """

    def build(self):
        # get devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")

        # arguments
        self.setattr_argument("channel_num",                        NumberValue(default=0, ndecimals=0, step=1, min=0, max=15))
        self.setattr_argument("voltage_v",                          NumberValue(default=5, ndecimals=3, step=1, min=0, max=10))
        self.setattr_argument("time_hold_ms",                       NumberValue(default=100, ndecimals=3, step=1, min=0.001, max=10000000000))

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.fastino0.set_dac(self.channel, self.voltage)
        delay(self.time_hold_ms * ms)
        self.fastino0.set_dac(self.channel, 0)

    def analyze(self):
        pass
