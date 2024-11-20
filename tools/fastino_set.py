from artiq.experiment import *


class FastinoSet(EnvExperiment):
    """
    Tool: Fastino Set

    Set a value on the Fastino.
    Useful for directly accessing the Fastino or when chaining a series of experiments.
    """

    def build(self):
        # get devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")

        # arguments
        self.setattr_argument("channel", NumberValue(default=0, precision=0, step=1, min=0, max=32))
        self.setattr_argument("voltage", NumberValue(default=0, precision=3, step=1, min=-10, max=10))

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.fastino0.set_dac(self.channel, self.voltage)

    def analyze(self):
        pass
