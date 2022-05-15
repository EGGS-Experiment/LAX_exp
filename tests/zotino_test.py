from artiq.experiment import *
import numpy as np
import time

class ZotinoTest(EnvExperiment):
    """Zotino Test"""

    def build(self):
        #self.setattr_argument("num_samples", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("delay_time", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("record_channel", NumberValue(ndecimals=0, step=1))

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
        self.urukul2_ch0.set_att(0)

        #set frequency and amplitude
        self.urukul2_ch0.set(1e8)
        self.urukul2_ch0.set_asf(0x3fff)

