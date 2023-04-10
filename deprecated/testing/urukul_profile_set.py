import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_DEST_FTW, RAM_MODE_DIRECTSWITCH


class urukul_profile_switch(EnvExperiment):
    """
    Urukul Profile Switching
    Switch DDS waveform quickly by using profiles
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        self.setattr_device("ttl8")

    def prepare(self):
        self.dds = self.urukul0_ch3

        # waveform
        self.waveform_asf = self.dds.amplitude_to_asf(0.7)

        # frequencies
        self.waveform_freq0 = self.dds.frequency_to_ftw(100 * MHz)
        self.waveform_freq1 = self.dds.frequency_to_ftw(200 * MHz)

        # timing
        self.time_delay_mu = self.core.seconds_to_mu(500 * ms)
        self.time_profile_mu = self.core.seconds_to_mu(1 * ms)
        self.time_ram_mu = self.core.seconds_to_mu(10 * ms)

        self.asf = self.urukul0_ch0.amplitude_to_asf(0.8)

        self.ftw = self.urukul0_ch0.frequency_to_ftw(50 * MHz)
        self.ftw2 = self.urukul0_ch0.frequency_to_ftw(60 * MHz)
        self.ftw3 = self.urukul0_ch0.frequency_to_ftw(55 * MHz)
        self.ftw4 = self.urukul0_ch0.frequency_to_ftw(100 * MHz)

    @kernel
    def run(self):
        self.core.reset()

        self.core.break_realtime()

        # self.urukul0_ch2.sw.on()
        # self.urukul0_ch3.sw.on()
        # self.core.break_realtime()
        # self.urukul0_ch2.set_mu(self.ftw, asf=self.asf, profile=0)
        # self.core.break_realtime()
        # self.urukul0_ch3.set_mu(self.ftw2, asf=self.asf, profile=0)
        #
        # delay(1*s)
        #
        # self.urukul0_ch2.set_mu(self.ftw3, asf=self.asf, profile=1)
        # self.core.break_realtime()
        #self.urukul0_ch3.set_mu(self.ftw4, asf=self.asf, profile=1)
        #self.ttl8.off()
        with parallel:
            self.urukul0_cpld.set_profile(1)
            self.ttl8.on()

        # self.dds.cpld.set_profile(0)
        # #self.dds.set_mu(self.ftw, asf=0x3fff, profile=0)
        # delay(3 * s)
        #
        #
        # self.dds.cpld.set_profile(1)
        # #self.dds.set_mu(self.ftw2, asf=0x1fff, profile=1)
        # #self.dds.sw.on()
        # self.dds.set_mu(self.ftw3, asf=self.asf, profile=2)
        # self.dds.cpld.set_profile(2)
        # delay(3 * s)
        #self.dds.set_mu(self.ftw4, asf=self.asf, profile=3)
        # self.dds.cpld.set_profile(3)
        # delay(3*s)


    def analyze(self):
        pass
