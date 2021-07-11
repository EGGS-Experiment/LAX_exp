from artiq.experiment import *

class MgmtTutorial(EnvExperiment):
    """Management Tutorial"""
    def build(self):
        pass #nothing
        #self.setattr_argument("count", NumberValue(ndecimals=0,step=1))

    def run(self):
        print("Hello World!")
        #for i in range(self.count):
        #   print("Hello World: ", i)