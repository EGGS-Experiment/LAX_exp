from artiq.experiment import *
from artiq.coredevice.ad9910 import _AD9910_REG_FTW

import numpy as np

class sampler_exp(EnvExperiment):
    """
    testing sampler
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device('ttl4')
        self.setattr_device('urukul0_ch0')
        self.setattr_device('urukul0_ch1')
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_ch0')
        self.setattr_device('urukul1_cpld')
        self.setattr_device('sampler0')

    @kernel
    def run(self):
        self.core.reset()
        self.sampler0.init()
        self.sampler0.set_gain_mu(0, 0)
        self.sampler0.set_gain_mu(1, 1)
        self.sampler0.set_gain_mu(2, 2)
        self.sampler0.set_gain_mu(3, 3)
        self.core.reset()
        print(self.sampler0.get_gains_mu())
        self.core.reset()
        yzde = self.sampler0.sample_mu([0]*8)
        self.core.reset()
        for val in yzde:
            print(val)

