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
        self.setattr_device("scheduler")
        self.setattr_device("ccb")

    def prepare(self):
        pass

    def run(self):
        print('\t\t\tLASERSCAN - RUNNING')
        for i in range(10):
            sleep(0.2)

    def analyze(self):
        res_dj = [[103.201, 0.5]]
        self.set_dataset('tmpres.ls.results', res_dj, broadcast=True, persist=False, archive=False)
        self.set_dataset('tmpres.ls.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
