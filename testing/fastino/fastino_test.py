from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class FastinoSet(EnvExperiment):
    """
    Fastino Set
    Set a value on the Fastino.
    """

    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # get dynamic devices
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul0_ch3')
        self.setattr_device("phaser0")
        self.setattr_device("fastino0")

        # get triggering TTLs
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

        # arguments
        self.setattr_argument("channel", NumberValue(default=0, precision=0, step=1, min=0, max=32))
        self.setattr_argument("voltage", NumberValue(default=0, precision=3, step=1, min=-10, max=10))

    def prepare(self):
        pass

    @kernel
    def run(self):
        # reset
        self.core.reset()
        self.core.break_realtime()

        # run
        for i in range(100000):
            with parallel:
                with sequential:
                    self.ttl8.on()
                    delay_mu(1830)
                    self.ttl9.on()
                self.fastino0.set_dac(0, 1.)

            delay_mu(100000)

            with parallel:
                self.ttl8.off()
                self.ttl9.off()

            self.fastino0.set_dac(0, 0.)
            delay_mu(100000)


        # clean up
        self.core.break_realtime()
        with parallel:
            self.ttl8.off()
            self.ttl9.off()
            self.fastino0.set_dac(0, 0.)


    def analyze(self):
        pass
