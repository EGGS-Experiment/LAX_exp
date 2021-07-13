from artiq.experiment import *
import numpy as np
import time

class UrukulTest(EnvExperiment):
    """Management Tutorial"""
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl4")
        self.setattr_device("urukul2_ch1")


    def run(self):
        self.core.reset()
        self.urukul2_ch1.cfg_sw(0, 1)
        self.urukul2_ch1.init()