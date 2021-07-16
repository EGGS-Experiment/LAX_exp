from artiq.experiment import *
import numpy as np
import time

class UrukulTest(EnvExperiment):
    """Urukul Test"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl4")
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")

    @kernel
    def run(self):
        self.core.reset()
        #initialize urukul
        self.urukul2_cpld.init()

        #initialize channel
        self.urukul2_ch0.cfg_sw(1)
        self.urukul2_ch0.init()

        #set frequency and amplitude
        self.urukul2_ch0.set(1e8)

