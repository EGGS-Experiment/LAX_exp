import numpy as np
from artiq.experiment import *

from time import sleep
from numpy import random


class autocalib_sbc_test(EnvExperiment):
    """
    autocalib sbc test
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
        # print('\t\t\tSBC - Test')
        res_dj = [[0.001, 0.12], [[1,2,3], [10,20,30]], [[4,5,6], [40,50,60]]]
        self.set_dataset('tmpres.sbc.result', res_dj, broadcast=True, persist=False, archive=False)
        self.set_dataset('tmpres.sbc.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
