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
        # self.setattr_device("core_dma")
        # self.setattr_device("ttl8")
        # self.setattr_device("ttl9")
        #
        # self.setattr_device("ttl0_counter")
        # self.setattr_device('urukul2_ch2')
        # self.setattr_device('urukul2_cpld')
        self.setattr_device('phaser0')

        # self.setattr_argument("thkim2", Scannable(
        #     default=[
        #         CenterScan(777.5, 4, 0.1, randomize=True),
        #         ExplicitScan([777.5]),
        #         ExplicitScan([767.2, 319.2, 1582, 3182]),
        #     ],
        #     global_min=0, global_max=10000, global_step=1,
        #     unit="kHz", scale=1, precision=3
        # ), group='EGGS_Heating.frequencies')

        # print(self.get_device_db())

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
        pass
        # self.cxn=labrad.connect()
        # self.dc=self.cxn.dc_server
        # self.voltage_set(24, 0.)
        # self.dc.toggle(24,1)


    @kernel
    def run(self):
        self.core.break_realtime()
        self.phaser0.init()
        self.core.break_realtime()

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att_mu(0xFF)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_att_mu(0xFF)
        self.core.break_realtime()

    # @kernel
    # def run_tmp(self):
    #     self.core.reset()
    #     self.core.break_realtime()
    #
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_nco_frequency(100. * MHz)
    #     # delay_mu(40)
    #     # self.phaser0.channel[1].set_nco_frequency(100. * MHz)
    #     # # clear DAC NCO phase
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_nco_phase(0.)
    #     # delay_mu(40)
    #     # self.phaser0.channel[1].set_nco_phase(0.)
    #     # # sync DAC for both channels
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.dac_sync()
    #     # # add slack
    #     # self.core.break_realtime()
    #
    #
    #
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].set_att(0. * dB)
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[1].set_att(0. * dB)
    #
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_duc_frequency(5. * MHz)
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[1].set_duc_frequency(5. * MHz)
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_duc_cfg()
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[1].set_duc_cfg()
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.duc_stb()
    #
    #
    #
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].oscillator[0].set_frequency(5. * MHz)
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0., clr=0)
    #     self.core.break_realtime()
    #
    #
    #     # self.core.reset()
    #     # self.ttl8.off()
    #     # self.ttl9.off()
    #     # self.core.break_realtime()
    #     # delay_mu(10000000)
    #     # self.core.break_realtime()
    #
    #     # self.ttl8.on()
    #     # self.core.wait_until_mu(now_mu())
    #     # delay_mu(500000)
    #     # self.voltage_set(24, 1.0)
    #     # self.ttl9.on()
    #     # self.core.wait_until_mu()
    #
    #
    #     # self.core.wait_until_mu(now_mu())
    #     # with parallel:
    #     # self.ttl8.on()
    #     # self.voltage_set(24, 1.0)
    #     # delay_mu(1000000)
    #     # self.core.wait_until_mu(now_mu())
    #     # delay_mu(10000000)
    #     # self.ttl8.on()
    #     # self.core.break_realtime()
    #
    #     # delay_mu(10000000)
    #     # self.ttl9.on()
    #
    #     # delay_mu(500000000)
    #     # self.ttl8.off()
    #     # self.ttl9.off()
