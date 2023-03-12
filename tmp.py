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

        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("ttl0_counter")

        self.repetitions = 200

        # self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        # self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
        #                                                                     default=CenterScan(85, 5, 0.2, randomize=True),
        #                                                                     global_min=30, global_max=400, global_step=1,
        #                                                                     unit="MHz", scale=1, ndecimals=5
        #                                                                 ))

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

        calib_timestamp = datetime.timestamp(datetime.now())
        # th0 = np.arange(85,137,2)
        # th1 = np.array([0.15625, 0.15625, 0.140625, 0.125, 0.1171875, 0.109375, 0.109375,
        #                 0.109375, 0.1171875, 0.1171875, 0.109375, 0.109375, 0.109375,
        #                 0.1171875, 0.125, 0.125, 0.125, 0.1328125, 0.140625, 0.140625,
        #                 0.15625, 0.171875, 0.203125, 0.25, 0.28125, 0.34375]) * 100
        #
        th0 = np.linspace(85,140,56)
        th1 = np.array([0.40625, 0.34375, 0.296875, 0.265625, 0.234375, 0.21875, 0.203125, 0.1875, 0.1796875, 0.171875, 0.171875, 0.1640625, 0.16015625, 0.15625, 0.15625, 0.1484375, 0.1484375, 0.1484375, 0.140625, 0.140625, 0.140625, 0.140625, 0.140625, 0.140625, 0.140625, 0.14453125, 0.14453125, 0.1484375, 0.1484375, 0.1484375, 0.15234375, 0.15625, 0.15625, 0.15625, 0.1640625, 0.1640625, 0.1640625, 0.171875, 0.171875, 0.171875, 0.1796875, 0.1875, 0.1875, 0.1953125, 0.203125, 0.21875, 0.2265625, 0.2421875, 0.2578125, 0.2734375, 0.296875, 0.3125, 0.34375, 0.375, 0.390625, 0.421875])*100
        # print(np.array([th0,th1]))
        self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', np.array([th0, th1]).transpose(), broadcast=True, persist=True)
        self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

    def prepare(self):
        pass
        # self.set_dataset('results', np.zeros((self.repetitions, 3)), broadcast=False)
        # self.set_dataset('results2', np.zeros((self.repetitions, 3)), broadcast=False)
        # self._iter_dataset = 0
        # self._iter_dataset2 = 0
        pass
        # self.setattr_device('urukul0_cpld')
        # self.setattr_device('urukul0_ch0')
        # self.delayth = np.int64(100)


    #@kernel(flags={"fast-math"})
    def run(self):
        pass
        # self.core.reset()
        # self.recordDMA()
        # self.core.break_realtime()
        # handle1 = self.core_dma.get_handle('tmphandle')
        # handle2 = self.core_dma.get_handle('tmphandle2')
        # self.core.break_realtime()
        # for i in range(self.repetitions):
        #     self.core_dma.playback_handle(handle1)
        #     self.core.break_realtime()
        #     self.core_dma.playback_handle(handle2)
        #     counts1 = self.ttl0_counter.fetch_count()
        #     self.core.break_realtime()
        #     counts2 = self.ttl0_counter.fetch_count()
        #
        #
        #     with parallel:
        #         with sequential:
        #             self.update_dataset(counts1)
        #             self.update_dataset2(counts2)
        #         self.core.break_realtime()

    @kernel
    def recordDMA(self):
        with self.core_dma.record('tmphandle'):

            self.urukul1_cpld.set_profile(1)

            # doppler cooling
            self.urukul1_cpld.cfg_switches(0b0110)
            delay_mu(3000000)
            self.urukul1_cpld.cfg_switches(0b0100)

            # switch to probe waveform
            self.urukul1_cpld.set_profile(0)

            # record fluorescence of probe beam
            self.urukul1_cpld.cfg_switches(0b0110)
            self.ttl0_counter.gate_rising_mu(9000000)
            self.urukul1_cpld.cfg_switches(0b0100)

        with self.core_dma.record('tmphandle2'):

            self.urukul1_cpld.set_profile(1)


            # doppler cooling
            self.urukul1_cpld.cfg_switches(0b0110)
            delay_mu(3000000)
            self.urukul1_cpld.cfg_switches(0b0100)

            # switch to probe waveform
            self.urukul1_cpld.set_profile(1)

            # record fluorescence of probe beam
            self.urukul1_cpld.cfg_switches(0b0000)
            self.ttl0_counter.gate_rising_mu(1000)
            self.urukul1_cpld.cfg_switches(0b0000)

    @rpc(flags={"async"})
    def print_rpc(self, text):
        print(text)

    @rpc(flags={"async"})
    def update_dataset(self, counts):
        self.mutate_dataset('results', self._iter_dataset, np.array([0, 0, counts]))
        self._iter_dataset += 1

    @rpc(flags={"async"})
    def update_dataset2(self, counts):
        self.mutate_dataset('results2', self._iter_dataset2, np.array([0, 0, counts]))
        self._iter_dataset2 += 1

    def analyze(self):
        print('test done')
