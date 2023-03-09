import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
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
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        self.setattr_device("ttl10")    # rf blank
        self.setattr_device("ttl11")    # oscilloscope trigger
        self.setattr_device("ttl12")    # rf switch

        # self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        # self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
        #                                                                     default=CenterScan(85, 5, 0.2, randomize=True),
        #                                                                     global_min=30, global_max=400, global_step=1,
        #                                                                     unit="MHz", scale=1, ndecimals=5
        #                                                                 ))

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

        calib_timestamp = datetime.timestamp(datetime.now())
        th0 = np.arange(85,135,2)
        th1 = np.array([0.15625, 0.140625, 0.125, 0.1171875, 0.109375, 0.109375,
                        0.109375, 0.1171875, 0.1171875, 0.109375, 0.109375, 0.109375,
                        0.1171875, 0.125, 0.125, 0.125, 0.1328125, 0.140625, 0.140625,
                        0.15625, 0.171875, 0.203125, 0.25, 0.28125, 0.34375]) * 100

        self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', np.array([th0, th1]).transpose(), broadcast=True, persist=True)
        self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

    def prepare(self):
        pass
        # self.setattr_device('urukul0_cpld')
        # self.setattr_device('urukul0_ch0')
        # self.delayth = np.int64(100)


    # @kernel(flags={"fast-math"})
   # @kernel
    def run(self):
        pass
        # self.core.reset()
        # #self.core.break_realtime()
        #
        # #at_mu(now_mu())
        # self.urukul0_ch0.set_mu(472446402, asf=3276)
        # #delay_mu(10)
        # self.urukul0_ch0.set_mu(515396075, asf=3276)
        # #delay_mu(10)
        # self.urukul0_ch0.set_mu(558345748, asf=3276)
        # #delay_mu(10)
        #
        # self.core.break_realtime()

    @kernel
    def recordDMA(self):
        pass

    @rpc(flags={"async"})
    def print_rpc(self, text):
        print(text)

    def analyze(self):
        print('test done')
