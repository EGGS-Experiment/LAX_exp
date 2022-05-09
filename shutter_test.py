from artiq.experiment import *

class ShutterTest(EnvExperiment):
    """Shutter Test"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl4")

    @kernel
    def run(self):
        self.core.reset()
        self.ttl4.output()
        for _ in 10:
            self.ttl4.on()
            delay(1)
            self.ttl4.off()
            delay(1)
