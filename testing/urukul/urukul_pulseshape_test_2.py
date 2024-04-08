from artiq.experiment import *  # Imports everything from experiment library
from artiq.coredevice.ad9910 import _AD9910_REG_CFR1, RAM_DEST_ASF, RAM_MODE_BIDIR_RAMP, RAM_MODE_CONT_RAMPUP, RAM_MODE_RAMPUP

import numpy as np


# This code demonstrates use of the urukul RAM. It produces a 125MHz pulse that ramps up in amplitude, holds a fixed amplitude and then ramps back down

class AD9910RAM2(EnvExperiment):
    '''Urukul RAM Amplitude Ramp 2'''

    def build(self):  # this code runs on the host computer
        self.setattr_device("core")  # sets core device drivers as attributes
        self.setattr_device("ttl8")  # sets ttl channel 6 device drivers as attributes
        self.setattr_device("ttl9")  # sets ttl channel 6 device drivers as attributes
        self.dds = self.get_device("urukul0_ch3")  # sets urukul 0, channel 1 device drivers as attributes and renames object self.dds


    def prepare(self):
        # self.amp = [0., 0.5, 0., 0.5, 0.]
        self.amp = [0.5, 0., 0., 0.5, 0., 0., 0.5]
        self.asf_ram = [0] * len(self.amp)

    # @kernel
    # def configure_ram_mode(self, dds):
    #     self.core.break_realtime()
    #
    #     dds.set_cfr1(ram_enable=0)
    #     self.cpld.io_update.pulse_mu(8)
    #     self.cpld.set_profile(0)  # Enable the corresponding RAM profile
    #     # Profile 0 is the default
    #     dds.set_profile_ram(start=0, end=len(self.asf_ram) - 1,
    #                         step=250, profile=0, mode=RAM_MODE_CONT_RAMPUP)
    #     self.cpld.io_update.pulse_mu(8)
    #     dds.amplitude_to_ram(self.amp, self.asf_ram)
    #     dds.write_ram(self.asf_ram)
    #     self.core.break_realtime()
    #     dds.set(frequency=5 * MHz, ram_destination=RAM_DEST_ASF)
    #     # Pass osk_enable=1 to set_cfr1() if it is not an amplitude RAM
    #     dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
    #     self.cpld.io_update.pulse_mu(8)

    @kernel
    def run(self):
        self.core.reset()

        self.dds.init()
        self.core.break_realtime()
        self.core.break_realtime()

        self.dds.cfg_sw(True)
        self.dds.set_att(8 * dB)
        self.core.break_realtime()
        delay_mu(1000000)
        at_mu(now_mu())

        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.cpld.set_profile(2)
        self.dds.set_profile_ram(start=0, end=len(self.asf_ram) - 1, step=2, profile=2, mode=RAM_MODE_RAMPUP)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.amplitude_to_ram(self.amp, self.asf_ram)
        self.dds.write_ram(self.asf_ram)
        self.core.break_realtime()

        self.dds.set(frequency=82 * MHz, ram_destination=RAM_DEST_ASF, profile=2)
        # self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF, phase_autoclear=1)
        self.dds.write32(_AD9910_REG_CFR1, (1 << 31) | (RAM_DEST_ASF << 29) | (1 << 16) | (1 << 13))

        with parallel:
            self.dds.cpld.io_update.pulse_mu(8)
            with sequential:
                delay_mu(80)
                self.ttl8.on()



        # cleanup
        at_mu(now_mu() + 1000000)
        self.ttl8.off()
        self.ttl9.off()
