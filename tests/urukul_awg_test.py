from artiq.experiment import *
from numpy import int32, int64


class Urukul_AWG(EnvExperiment):
    """
    Urukul AWG Test.
    """

    def build(self):
        """
        Ensure that the necessary devices are set.
        """
        self.setattr_device("core")                                 # always needed
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")

    def prepare(self):
        # set values
        # get amplitude RAM data
        pass

    @kernel
    def run(self):
        self.core.reset()
        # set amplitude awg mode
        # program RAM
        # set up config registers correctly
        # run pulse

    def analyze(self):
        pass