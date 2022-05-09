from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class DMARun(EnvExperiment):
    """
    Run a DMA Exp by handle.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        #ttl0 count read
        #ttl1 over light
        #ttl2 linetrigger in
        #ttl4 power
        #ttl5 linetrigger out

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.core_dma.playback_handle(self.core_dma.get_handle('PMT_exp'))

    def analyze(self):
        th1 = self.get_dataset('pmt_test_dataset')
        print(th1)
