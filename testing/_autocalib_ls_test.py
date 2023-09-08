import numpy as np
from artiq.experiment import *

from time import sleep
from numpy import random


class autocalib_ls_test(EnvExperiment):
    """
    autocalib ls test
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("scheduler")
        self.setattr_device("ccb")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

    def prepare(self):
        pass

    def run(self):
        print('idk')
        for i in range(10):
            sleep(0.5)

    def analyze(self):
        print('\t\t\tLS - TEST')
