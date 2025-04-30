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

    def prepare(self):
        self.dds = self.get_device('urukul0_ch1')
        self.ftw = self.dds.frequency_to_ftw(120.339 * MHz)
        self.asf = self.dds.amplitude_to_asf(0.5)
        self.att_db = 7 * dB

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

