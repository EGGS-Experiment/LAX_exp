from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class PIDTest(EnvExperiment):
    """
    Test PID.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        # set up
        self.
        # loop
        # read in pickoff voltage
        # read in setpoint voltage
        # do PID
        # output error signal voltage
        self.fastino0.set_dac(0, 0)

    def analyze(self):
        pass
