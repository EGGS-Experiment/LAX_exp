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
        self.setattr_device('urukul1_ch2')
        self.setattr_device('urukul1_cpld')

        # urukul devices
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

        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")
        self.setattr_device("urukul2_ch2")
        self.setattr_device("urukul2_ch3")

    def prepare(self):
        self.dds =      self.urukul1_ch2
        self.dds_cpld = self.dds.cpld

        self.freq =     np.int32(0)
        self.ampl =     np.int32(0)
        self.phase =    np.int32(0)
        self.att =      np.int32(0)

        self._profile0_word = np.int64(0)
        self._profile1_word = np.int64(0)
        self._profile2_word = np.int64(0)
        self._profile3_word = np.int64(0)
        self._profile4_word = np.int64(0)
        self._profile5_word = np.int64(0)
        self._profile6_word = np.int64(0)
        self._profile7_word = np.int64(0)

        self.th0 = (0., 0., 0.)
        self.th1 = (0., 0., 0.)


    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()

        # tmp remove
        self.dds_cpld.init()
        self.core.break_realtime()
        self.core.break_realtime()
        delay_mu(10000000)
        self.core.break_realtime()
        # tmp remove

        # tmp remove
        self.dds.init()
        self.core.break_realtime()
        self.core.break_realtime()
        delay_mu(10000000)
        self.core.break_realtime()
        # tmp remove

        # tmp remove
        self.dds.set_mu(0x1F0FFFF, asf=0xFFF, profile=1)
        self.core.break_realtime()
        # tmp remove


        self.dds_cpld.set_profile(0)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self._profile0_word = np.int64(self.dds.read64(0x0E))
        self.core.break_realtime()

        self.dds_cpld.set_profile(1)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self._profile1_word = np.int64(self.dds.read64(0x0F))
        self.core.break_realtime()

        self.dds_cpld.set_profile(2)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self._profile2_word = np.int64(self.dds.read64(0x10))
        self.core.break_realtime()

        self.dds_cpld.set_profile(3)
        self.core.break_realtime()
        self._profile3_word = np.int64(self.dds.read64(0x11))
        self.core.break_realtime()

        self.dds_cpld.set_profile(4)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self._profile4_word = np.int64(self.dds.read64(0x12))
        self.core.break_realtime()

        self.dds_cpld.set_profile(7)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()
        self._profile7_word = np.int64(self.dds.read64(0x15))
        self.core.break_realtime()

        self.th0 = self.dds.get(profile=0)
        self.core.break_realtime()

        self.th1 = self.dds.get(profile=2)
        self.core.break_realtime()


    def analyze(self):
        ftw_0 = np.int32(self._profile0_word & 0xFFFFFFFF)
        pow_0 = np.int32((self._profile0_word >> 32) & 0xFFFF)
        asf_0 = np.int32((self._profile0_word >> 48) & 0x3FFF)

        ftw_1 = np.int32(self._profile1_word & 0xFFFFFFFF)
        pow_1 = np.int32((self._profile1_word >> 32) & 0xFFFF)
        asf_1 = np.int32((self._profile1_word >> 48) & 0x3FFF)

        ftw_2 = np.int32(self._profile2_word & 0xFFFFFFFF)
        pow_2 = np.int32((self._profile2_word >> 32) & 0xFFFF)
        asf_2 = np.int32((self._profile2_word >> 48) & 0x3FFF)

        ftw_3 = np.int32(self._profile3_word & 0xFFFFFFFF)
        pow_3 = np.int32((self._profile3_word >> 32) & 0xFFFF)
        asf_3 = np.int32((self._profile3_word >> 48) & 0x3FFF)

        ftw_4 = np.int32(self._profile4_word & 0xFFFFFFFF)
        pow_4 = np.int32((self._profile4_word >> 32) & 0xFFFF)
        asf_4 = np.int32((self._profile4_word >> 48) & 0x3FFF)

        ftw_7 = np.int32(self._profile7_word & 0xFFFFFFFF)
        pow_7 = np.int32((self._profile7_word >> 32) & 0xFFFF)
        asf_7 = np.int32((self._profile7_word >> 48) & 0x3FFF)

        print(self.th0)
        print(self.th1)



        print("\tprofile 0:")
        print("\t\tfreq: {:.4f} MHz".format(self.dds.ftw_to_frequency(ftw_0) / MHz))
        print("\t\tampl: {:.4f} %".format(self.dds.asf_to_amplitude(asf_0) * 100.))
        print("\t\tphas: {:.4f} turns\n".format(self.dds.pow_to_turns(pow_0)))

        print("\tprofile 1:")
        print("\t\tfreq: {:.4f} MHz".format(self.dds.ftw_to_frequency(ftw_1) / MHz))
        print("\t\tampl: {:.4f} %".format(self.dds.asf_to_amplitude(asf_1) * 100.))
        print("\t\tphas: {:.4f} turns\n".format(self.dds.pow_to_turns(pow_1)))

        print("\tprofile 2:")
        print("\t\tfreq: {:.4f} MHz".format(self.dds.ftw_to_frequency(ftw_2) / MHz))
        print("\t\tampl: {:.4f} %".format(self.dds.asf_to_amplitude(asf_2) * 100.))
        print("\t\tphas: {:.4f} turns\n".format(self.dds.pow_to_turns(pow_2)))

        print("\tprofile 3:")
        print("\t\tfreq: {:.4f} MHz".format(self.dds.ftw_to_frequency(ftw_3) / MHz))
        print("\t\tampl: {:.4f} %".format(self.dds.asf_to_amplitude(asf_3) * 100.))
        print("\t\tphas: {:.4f} turns\n".format(self.dds.pow_to_turns(pow_3)))

        print("\tprofile 7:")
        print("\t\tfreq: {:.4f} MHz".format(self.dds.ftw_to_frequency(ftw_7) / MHz))
        print("\t\tampl: {:.4f} %".format(self.dds.asf_to_amplitude(asf_7) * 100.))
        print("\t\tphas: {:.4f} turns\n".format(self.dds.pow_to_turns(pow_7)))

