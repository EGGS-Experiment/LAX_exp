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





        # self.urukul2_ch2.cfg_sw(False)
        # self.urukul2_ch2.sw.off()
        # self.core.break_realtime()
    def prepare(self):
        _CONFIG = "th0_terminated"

        import labrad
        self.cxn=labrad.connect()
        self.na=self.cxn.network_analyzer_server
        self.na.select_device()
        self.num_points = self.na.sweep_points(1601)
        self.freq_start = self.na.frequency_start()
        self.freq_stop = self.na.frequency_stop()
        self._results0 = np.linspace(self.freq_start, self.freq_stop, self.num_points)
        self._results1 = np.zeros(self.num_points)

        self.na.gpib_write('DISP:WIND1:TRAC:Y:AUTO ONCE')

        self.set_dataset('config', _CONFIG, broadcast=False)
        self.set_dataset('freq_start_hz', self.freq_start, broadcast=False)
        self.set_dataset('freq_stop_hz', self.freq_stop, broadcast=False)
        self.set_dataset('num_points', self.num_points, broadcast=False)

        # self.set_dataset('results', np.zeros((self.repetitions, 3)), broadcast=False)
        # self.set_dataset('results2', np.zeros((self.repetitions, 3)), broadcast=False)
        # self._iter_dataset = 0
        # self._iter_dataset2 = 0
        # yz0=self.get_dataset('calibration.temperature.asf_calibration_curve_mhz_pct')

    def run(self):
        res0 = self.na.gpib_query('TRAC? CH1FDATA')
        res1 = np.array([float(strval) for strval in res0.split(',')])
        self._results1 = res1
    #
    #
    def analyze(self):
        res_fin = np.array([self._results0, self._results1]).transpose()
        self.set_dataset('results', res_fin)
        print(res_fin)
    #
