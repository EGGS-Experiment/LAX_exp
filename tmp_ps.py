from artiq.experiment import *  # Imports everything from experiment library
from artiq.coredevice.ad9910 import (  # Imports RAM destination amplitude scale factor and RAM mode bidirectional ramp methods from AD9910 Source
    RAM_DEST_ASF, RAM_MODE_BIDIR_RAMP)

import numpy as np


# This code demonstrates use of the urukul RAM. It produces a 125MHz pulse that ramps up in amplitude, holds a fixed amplitude and then ramps back down

class AD9910RAM(EnvExperiment):
    '''Urukul RAM Amplitude Ramp'''

    def build(self):  # this code runs on the host computer
        self.setattr_device("core")  # sets core device drivers as attributes
        self.setattr_device("ttl15")  # sets ttl channel 6 device drivers as attributes
        self.dds = self.get_device("urukul0_ch1")  # sets urukul 0, channel 1 device drivers as attributes and renames object self.dds


    def prepare(self):
        # produce data to be loaded to RAM
        # self.data_len = 512
        # self.data = np.zeros(self.data_len, dtype=np.int32)
        # for i in range(self.data_len):
        #     ampl_pct = np.power(np.sin(i / self.data_len * np.pi), 2) / 2.1
        #     self.data[i] = self.dds.amplitude_to_asf(ampl_pct)
        #
        # self.data = list(self.data)

        #print(self.data)
        #raise Exception('thkim')
        # tmp 1
        self.freq_ftw = self.dds.frequency_to_ftw(103.771 * MHz)
        self.freq_2 = 51 * MHz


    @kernel  # this code runs on the FPGA
    def run(self):
        # data_len = 512
        # data = [0] * 512
        # for i in range(data_len):
        #     # ampl_pct = np.power(np.sin(i / self.data_len * np.pi), 2) / 2.1
        #     # self.data[i] = self.dds.amplitude_to_asf(0.5)
        #     data[i] = 0xfff

        n = 8                                                              #defines variable n for list length exponent
        data = [0]*(1 << n)                                                 #declares list as 2^n integer values
        for i in range(len(data)//2):                                       #splits list into 2 halves and defines each separately
            data[i] = i << (32 - (n - 1))                                   #first half ramps up to maximum amplitude in machine units
            data[i + len(data)//2] = 0xfff << 16                           #second half holds maximum amplitude

        # reset core
        self.core.reset()
        self.dds.init()
        self.core.break_realtime()
        self.core.break_realtime()
        self.dds.cpld.init()
        self.core.break_realtime()
        self.core.break_realtime()

        # set ram profile
        self.dds.set_profile_ram(  # sets profile in RAM to be used
            start=0, end=0+len(data) - 1, step=10,
            # start/end give addresses of ends of ram data, step gives step length
            profile=0, mode=RAM_MODE_BIDIR_RAMP)  # mode: bidirectional ramp

        self.dds.cpld.set_profile(0)  # sets CPLD profile pins
        self.dds.cpld.io_update.pulse_mu(8)  # I think this clocks all the CPLD registers so they take the values written to them
        delay(1 * ms)  # 1ms delay

        # write to ram
        self.dds.write_ram(data)  # writes data list to ram
        delay(10 * ms)  # 10ms delay

        # write to cfr
        self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)  # writes to CFR1 (control function register 1)
        # enables ram, sets ram data as amplitude scale factor

        # set urukuln parameters and turn channel on
        self.dds.set_frequency(self.freq_2)
        self.dds.cpld.io_update.pulse_mu(8)  # I think this clocks all the CPLD registers so they take the values written to them
        self.dds.set_att(10 * dB)
        self.dds.cfg_sw(True)
        self.core.break_realtime()

        # thkim
        self.dds.set_mu(self.freq_ftw, asf=0x1fff, profile=7)
        self.core.break_realtime()


        for i in range(1000):  # loops until manually broken
            self.core.break_realtime()
            self.dds.cfg_sw(True)

            self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
            self.dds.cpld.io_update.pulse_mu(8)

            with parallel:  # runs indented code in parallel
                self.ttl15.pulse(1 * us)  # 1us TTL pulse for triggering oscilloscope
                self.dds.cpld.set_profile(0)  # profile 0 tells CPLD to start ramping up

            delay(2 * us)  # 2us delay

            with parallel:  # runs indented code in parallel
                self.ttl15.pulse(1 * us)  # 1us TTL pulse
                self.dds.cpld.set_profile(1)  # profile 1 tells CPLD to start ramping back down

            delay(1 * ms)  # 1ms delay
            self.dds.cfg_sw(False)
            self.core.break_realtime()

            self.dds.set_cfr1(ram_enable=0)
            self.dds.cpld.io_update.pulse_mu(8)
            self.core.break_realtime()
            #
            self.dds.cpld.set_profile(7)
            self.dds.cpld.io_update.pulse_mu(8)
            self.dds.cfg_sw(True)
            delay_mu(20000)
            self.dds.cfg_sw(False)
