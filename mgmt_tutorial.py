from artiq.experiment import *
import numpy as np
import time

class MgmtTutorial(EnvExperiment):
    """Management Tutorial"""
    def build(self):
        pass #nothing
        #self.setattr_argument("count", NumberValue(ndecimals=0,step=1))


    def run(self):
        print("Hello World!")
        # self.set_dataset("parabola", np.full(self.count, np.nan), broadcast=True)
        # for i in range(self.count):
        #     self.mutate_dataset("parabola", i, i * i)