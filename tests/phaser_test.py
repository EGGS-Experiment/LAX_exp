from artiq.experiment import *
from numpy import int32, int64


class PhaserTest(EnvExperiment):
    """
    Phaser Test.
    Try to generate an 80MHz wave with an erf window just like we simulated in EGGS.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")

        self.setattr_device("fastino0")
        self.th1 = self.urukul0_ch0.frequency_to_ftw(100*MHz)
        self.th2 = self.urukul0_ch0.amplitude_to_asf(0.4)

    @kernel
    def run(self):
        self.core.reset()

    def analyze(self):
        pass
