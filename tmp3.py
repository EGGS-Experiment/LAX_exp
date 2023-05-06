import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testarg34(EnvExperiment):
    """
    testarg34
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        self.setattr_device("ttl10")    # rf blank
        self.setattr_device("ttl11")    # oscilloscope trigger
        self.setattr_device("ttl12")    # rf switch

        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("ttl0_counter")

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)

        # calib_timestamp = datetime.timestamp(datetime.now())
        # th0 = np.arange(85,137,2)
        # th1 = np.array([0.15625, 0.15625, 0.140625, 0.125, 0.1171875, 0.109375, 0.109375,
        #                 0.109375, 0.1171875, 0.1171875, 0.109375, 0.109375, 0.109375,
        #                 0.1171875, 0.125, 0.125, 0.125, 0.1328125, 0.140625, 0.140625,
        #                 0.15625, 0.171875, 0.203125, 0.25, 0.28125, 0.34375]) * 100

        # # print(np.array([th0,th1]))
        # self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', np.array([th0, th1]).transpose(), broadcast=True, persist=True)
        # self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

    def prepare(self):
        yz0=self.get_dataset('calibration.eggs.resonance_ratio_curve_mhz')
        print(yz0)

    def run(self):
        pass


    def analyze(self):
        pass
        # res_fin = np.array([self._results0, self._results1]).transpose()
        # self.set_dataset('results', res_fin)
        # print(res_fin)

