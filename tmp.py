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
        self.set_dataset('ampl_pump_readout_pct', 38.0, broadcast=True, persist=True)
        self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)
        self.set_dataset('ampl_repump_qubit_pct', 14.0, broadcast=True, persist=True)

        self.set_dataset('freq_pump_cooling_mhz', 110.0, broadcast=True, persist=True)
        self.set_dataset('freq_pump_readout_mhz', 110.0, broadcast=True, persist=True)
        self.set_dataset('freq_repump_cooling_mhz', 110.0, broadcast=True, persist=True)
        self.set_dataset('freq_repump_qubit_mhz', 110.0, broadcast=True, persist=True)

    def run(self):
        pass
