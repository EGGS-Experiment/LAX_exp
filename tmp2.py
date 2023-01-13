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
        # self.setattr_argument("test1", PYONValue({'a':1,'b':2}))
        # self.setattr_device('urukul1_ch0')

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    #def prepare(self):
        # self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        # self.os = self.cxn.oscilloscope_server
        # print(self.os.list_devices())

    #@kernel
    def run(self):
        self.set_dataset('beams.freq_mhz.freq_pump_rescue_mhz', 95.0, broadcast=True, persist=True)
        self.set_dataset('beams.ampl_pct.ampl_pump_rescue_pct', 50.0, broadcast=True, persist=True)
        # self.core.reset()
        # self.ttl8.on()
        # self.core.break_realtime()
        #pass
        # self.core.reset()
        # self.core.break_realtime()
        # self.urukul1_ch0.init()
        # delay(10*ms)
        #print('scde')


        #self.set_dataset('dds.dds_board_tickle_num', 0, broadcast=True, persist=True)
        #self.set_dataset('dds.dds_tickle_channel', 3, broadcast=True, persist=True)
        #self.set_dataset('pmt_gating_edge', 'rising', broadcast=True, persist=True)
        #print('type: {}'.format(type(self.test1)))
        #print(self.test1)

    def analyze(self):
        print('test done')
