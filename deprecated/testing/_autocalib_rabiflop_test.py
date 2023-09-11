import numpy as np
from artiq.experiment import *

from time import sleep
from numpy import random


class autocalib_rabiflop_test(EnvExperiment):
    """
    autocalib rabiflop test
    Testing
    """

    def build(self):
        self.setattr_device("scheduler")
        self.setattr_device("ccb")

    def prepare(self):
        pass

    def run(self):
        print('\t\t\tRABIFLOP - RUNNING')
        for i in range(10):
            sleep(0.2)

    def analyze(self):
        res_dj = [12.5, 0.1]
        self.set_dataset('tmpres.rabiflop.results', res_dj, broadcast=True, persist=False, archive=False)
        self.set_dataset('tmpres.rabiflop.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
