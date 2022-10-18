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

    def prepare(self):
        self.set_dataset('ampl_pump_cooling_pct', 22.8, broadcast=True, persist=True)

    def run(self):
        pass
