from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class FastinoSet(EnvExperiment):
    """
    Set a value on the Fastino.
    """

    def build(self):
        # get devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")

        # arguments
        self.setattr_argument("channel", NumberValue(default=0, ndecimals=0, step=1, min=0, max=32))
        self.setattr_argument("voltage", NumberValue(default=0, ndecimals=3, step=1, min=-10, max=10))

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.fastino0.set_dac(self.channel, self.voltage)

    def analyze(self):
        pass
