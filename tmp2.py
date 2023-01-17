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

    def prepare(self):
        self.hasenv = testhasenv(self, {'aa':11,'bbb':22}, 9031)
        # self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        # self.os = self.cxn.oscilloscope_server
        # print(self.os.list_devices())

    #@kernel
    def run(self):
        print('running')

    def analyze(self):
        print('test done')


class testhasenv(HasEnvironment):
    """
    test hasenv
    yzde
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_argument("test1", PYONValue({'a':1,'b':2}))
        self.setattr_argument("test2", NumberValue(default=104.335, ndecimals=5, step=1, min=1, max=10000))

    def prepare(self):
        print('\thasenv test:')
        print('\t\ttest1: {}'.format(self.test1))
        print('\t\ttest2: {}'.format(self.test2))
