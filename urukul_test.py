from artiq.experiment import *
import numpy as np
import time

class UrukulTest(EnvExperiment):
    """Urukul Test"""

    def build(self):
        #self.setattr_argument("frequency", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("amplitude", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("phase", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("channel", NumberValue(ndecimals=0, step=1))

        self.setattr_device("core")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")

    @kernel
    def run(self):
        self.core.reset()

        #initialize urukul
        self.urukul1_cpld.init(blind=True)
        self.urukul1_ch0.init()
        self.core.break_realtime()

        self.urukul1_ch0.cfg_sw(1)
        self.urukul1_ch0.set_att_mu(255)

        #set frequency and amplitude
        self.urukul1_ch0.set(1e8)
        self.urukul1_ch0.set_asf(0x3fff)