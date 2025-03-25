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
    kernel_invariants = {
        "dds"
    }

    def build(self):
        self.setattr_device("core")
        # self.setattr_device("ttl0_counter")
        # self.setattr_device('urukul2_ch2')
        # self.setattr_device('urukul2_cpld')


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
        self.dds = self.get_device('urukul0_ch1')
        self.ftw = self.dds.frequency_to_ftw(80 * MHz)
        self.asf = self.dds.amplitude_to_asf(0.58)
        self.att_db = 3 * dB

        self.profile_target = 7

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset
        self.core.reset()
        self.dds.cpld.io_update.pulse_mu(8)

        # keep board attenuations
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()

        # close switch while we update
        self.dds.sw.off()

        # configure
        self.dds.set_cfr1()
        self.dds.cpld.io_update.pulse_mu(8)

        # # set
        # self.dds.set_mu(self.ftw, asf=self.asf, profile=self.profile_target)
        # self.dds.cpld.io_update.pulse_mu(8)
        # self.dds.set_att(self.att_db * dB)
        # self.core.break_realtime()
        #
        # self.dds.cpld.set_profile(self.profile_target)
        # self.dds.cpld.io_update.pulse_mu(8)
        # self.core.break_realtime()
        #
        # # finish
        # self.dds.sw.on()

        # set ALL
        for i in range(8):
            self.core.break_realtime()
            self.dds.set_mu(self.ftw, asf=self.asf, profile=i)
            self.dds.cpld.io_update.pulse_mu(8)
            self.dds.set_att(self.att_db * dB)
            self.core.break_realtime()

            self.dds.cpld.set_profile(self.profile_target)
            self.dds.cpld.io_update.pulse_mu(8)
            self.core.break_realtime()

            # finish
            self.dds.sw.on()

