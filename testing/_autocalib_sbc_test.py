import numpy as np
from artiq.experiment import *

class autocalib_sbc_test(EnvExperiment):
    """
    autocalib sbc test
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

    def prepare(self):
        pass

    def run(self):
        print('\t\t\tSBCTEST')
