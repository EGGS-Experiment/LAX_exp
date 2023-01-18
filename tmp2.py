import numpy as np
import labrad
from os import environ
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
        self.setattr_device("ttl8")

        # self.setattr_device("core_dma")
        # self.setattr_device('urukul1_ch0')

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    def prepare(self):
        pass
        # self.call_child_method('prepare')
        # self.hasenv = testhasenv(self, {'aa':11,'bbb':22}, 9031)
        # self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        # self.os = self.cxn.oscilloscope_server
        # print(self.os.list_devices())

    #@kernel
    def run(self):
        self.set_dataset('pmt.input_channe', 'rising', broadcast=True, persist=True)

    def analyze(self):
        print('test done')
