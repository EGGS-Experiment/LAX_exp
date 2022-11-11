import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG


class Testing(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")

        self.time0 = self.core.seconds_to_mu(1000 * ms)
        self.time1 = self.core.seconds_to_mu(100 * us)

        self.set_dataset("tmp", [])
        self.setattr_dataset("tmp")

    @kernel
    def run(self):
        self.core.reset()

        self.ttl0.input()
        self.core.break_realtime()

        th1 = self.ttl0.count(self.ttl0.gate_rising_mu(self.time1))
        self.core.break_realtime()
        self.append_to_dataset("tmp", th1)
        #self.set_dataset('dds_tickle_channel', 3, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    def analyze(self):
        print(self.tmp)
