from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class FastinoSet(EnvExperiment):
    """Fastino Set Value"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.fastino0.set_dac(0, 1)

    def analyze(self):
        pass