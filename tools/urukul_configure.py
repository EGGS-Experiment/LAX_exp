import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class UrukulConfigure(EnvExperiment):
    """
    Tool: Urukul Configure

    Configure values for an urukul channel (AD9910 ONLY).
    """
    name = 'Urukul Configure'

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_target", EnumerationValue(list(dds_device_list), default='urukul2_ch1'))

        # CPLD parameters
        self.setattr_argument('dds_profile',        NumberValue(default=7, precision=0, step=1, min=0, max=7), group="dds.cpld")
        self.setattr_argument("switch_on",          BooleanValue(default=False), group="dds.cpld")
        self.setattr_argument('att_db',             NumberValue(default=14., precision=1, step=0.5, min=0., max=31.5), group="dds.cpld")

        # AD9910 parameters
        self.setattr_argument("initialize_dds", BooleanValue(default=False), group="dds.channel")
        self.setattr_argument('freq_mhz',       NumberValue(default=110, precision=7, step=5, min=1., max=400.), group="dds.channel")
        self.setattr_argument('ampl_pct',       NumberValue(default=50., precision=3, step=5., min=0., max=100.), group="dds.channel")
        self.setattr_argument('pow_turns',      NumberValue(default=0., precision=4, step=0.1, min=-1., max=1.), group="dds.channel")

        # todo: clear RAM mode via CFR1 (and maybe CFR2?)


    def prepare(self):
        self.dds = self.get_device('urukul0_ch1')
        self.ftw = self.dds.frequency_to_ftw(80 * MHz)
        self.asf = self.dds.amplitude_to_asf(0.61)
        self.att_db = 3. * dB

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

        # set
        self.dds.set_mu(self.ftw, asf=self.asf, profile=self.profile_target)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_att(self.att_db * dB)
        self.core.break_realtime()

        self.dds.cpld.set_profile(self.profile_target)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # finish
        self.dds.sw.on()

