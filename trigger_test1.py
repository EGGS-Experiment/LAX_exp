from artiq.experiment import *
from artiq.coredevice.ad9910 import _AD9910_REG_FTW

from numpy import zeros, ones
from datetime import datetime


class trigger_test(EnvExperiment):
    """
    testing linetrigger
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('ttl0')
        self.setattr_device('ttl1')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')
        self.set_dataset('ttt', zeros(1000), broadcast=True)
        self.setattr_dataset('ttt')
        self.i = 0


    @kernel
    def run(self):
        self.core.reset()
        self.ttl0.input()
        self.ttl1.input()
        self.core.reset()
        self.core.reset()
        handle = self.core_dma.get_handle('PMT_exp')
        self.core.reset()
        while True:
            if self.ttl0.count(self.ttl0.gate_falling(1 * ms)) > 0:
                break
            self.i += 1
            delay(4 * ms)
        self.core.reset()
        self.core_dma.playback_handle(handle)

    def analyze(self):
        print('total loops:', self.i)
