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

        self.setattr_device("ttl0_counter")
        self.setattr_device('urukul2_ch2')
        self.setattr_device('urukul2_cpld')

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
        self.cxn=labrad.connect()
        # self.dc=self.cxn.dc_server
        # self.voltage_set(24, 0.)
        # self.dc.toggle(24,1)

    @rpc
    def voltage_set(self, channel: TInt32, voltage_v: TFloat) -> TNone:
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        voltage_set_v = self.dc.voltage_fast(channel, voltage_v)
        # print('\tvoltage set: {}'.format(voltage_set_v))

    # @kernel
    def run(self):
        # self.core.reset()
        # self.ttl8.off()
        # self.ttl9.off()
        # self.core.break_realtime()
        # delay_mu(10000000)
        # self.core.break_realtime()

        # self.ttl8.on()
        # self.core.wait_until_mu(now_mu())
        # delay_mu(500000)
        # self.voltage_set(24, 1.0)
        # self.ttl9.on()
        # self.core.wait_until_mu()



        # self.core.wait_until_mu(now_mu())
        # with parallel:
        # self.ttl8.on()
        # self.voltage_set(24, 1.0)
        # delay_mu(1000000)
        # self.core.wait_until_mu(now_mu())
        # delay_mu(10000000)
        # self.ttl8.on()
        # self.core.break_realtime()

        # delay_mu(10000000)
        # self.ttl9.on()

        delay_mu(500000000)
        self.ttl8.off()
        self.ttl9.off()
