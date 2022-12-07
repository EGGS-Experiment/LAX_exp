import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testarg12(EnvExperiment):
    """
    testarg34
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_argument("test1", PYONValue([1,2,3]))

        # self.set_dataset('dds_tickle_channel', 3, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    #@kernel
    def run(self):
        #self.set_dataset('dds.dds_board_tickle_num', 0, broadcast=True, persist=True)
        #self.set_dataset('dds.dds_tickle_channel', 3, broadcast=True, persist=True)
        #self.set_dataset('pmt_gating_edge', 'rising', broadcast=True, persist=True)
        print(self.test1)
