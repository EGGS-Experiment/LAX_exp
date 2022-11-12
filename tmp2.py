import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

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

        self.time0 = self.core.seconds_to_mu(1 * ms)
        self.time1 = self.core.seconds_to_mu(10 * us)
        self.timed = self.core.seconds_to_mu(5 * us)

        self.set_dataset("tmp", [])
        self.setattr_dataset("tmp")

        # self.set_dataset('dds_tickle_channel', 3, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    #@kernel
    def run(self):
        print(self.core.seconds_to_mu(1 * us))
        #self.set_dataset('rf_mod_trigger_channel', 3, broadcast=True, persist=True)
