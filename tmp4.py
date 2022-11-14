import numpy as np
from artiq.experiment import *

class Testing(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl3")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.time0 = self.core.seconds_to_mu(1 * ms)
        self.time1 = self.core.seconds_to_mu(10 * us)
        self.timed = self.core.seconds_to_mu(5 * us)

        self.set_dataset("tmp", [])
        self.setattr_dataset("tmp")

    @kernel
    def run(self):
        self.core.reset()

        self.ttl0.input()
        self.ttl3.input()
        self.core.break_realtime()

        # record dma and get handle
        handle_initialize = self.DMArecord()
        self.core.break_realtime()

        for i in range(100):

            self.core_dma.playback_handle(handle_initialize)
            self.append_to_dataset("tmp", i)
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        with self.core_dma.record("yzde"):
            self.testfunc()
            delay_mu(self.timed)
            self.testfunc2()

        return self.core_dma.get_handle("yzde")

    @kernel
    def testfunc(self):
        self.urukul0_ch3.cfg_sw(0)
        delay_mu(self.timed)
        self.urukul0_ch3.cfg_sw(1)

    @kernel
    def testfunc2(self):
        self.urukul0_ch2.cfg_sw(0)
        delay_mu(self.timed)
        self.urukul0_ch2.cfg_sw(1)


    def analyze(self):
        print(self.tmp)
