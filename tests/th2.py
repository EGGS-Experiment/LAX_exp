from artiq.experiment import *
from artiq.coredevice.ad9910 import (RAM_DEST_ASF, RAM_DEST_FTW, RAM_MODE_BIDIR_RAMP)


class UrukulFrequeencyRamp(EnvExperiment):
    """
    Urukul Frequency Ramp
    """

    def build(self):
        self.setattr_device("core")
        self.u = self.get_device("urukul0_ch3")
        self.th = self.get_device('urukul0_ch0')

    @kernel
    def run(self):
        # produce data to be loaded to RAM
        n = 9
        data = [0] * (1 << n)
        # todo: 0x1c28f5c2 is 110mhz
        # todo: 0xccccccc is 50mhz
        for i in range(len(data) // 2):
            # first half ramps up to maximum amplitude in machine units
            data[i] = i << (32 - (n - 1))
            # second half holds maximum amplitude
            data[i + len(data) // 2] = 0xffff << 16

        self.core.reset()
        # initialize
        #self.dds.cpld.init()
        self.core.break_realtime()
        #self.dds.init()
        self.core.break_realtime()
        delay(1 * ms)

        # set ram profile
        self.dds.set_profile_ram(
            start=0, end=len(data) - 1, step=0x01,
            profile=0,
            mode=RAM_MODE_BIDIR_RAMP
        )

        # set CPLD profile pins
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay(10 * ms)

        # write to ram
        self.core.break_realtime()
        self.dds.write_ram(data)  # writes data list to ram
        self.core.break_realtime()
        self.core.break_realtime()
        self.core.break_realtime()

        # writes to CFR1 (control function register 1)
        # enables ram, sets ram data as amplitude scale factor
        self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_FTW)

        # set urukul parameters and turn channel on
        self.dds.set_frequency(50 * MHz)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_att(10 * dB)
        self.dds.sw.on()
        self.core.break_realtime()
        self.core.break_realtime()

        # tmp remove: set urukul0_ch0
        self.th.set_mu(472446402, asf=0x1fff)
        self.core.break_realtime()
        self.th.set_att(8 * dB)
        self.core.break_realtime()
        self.th.sw.on()
        self.core.break_realtime()

        # continually ramp
        while True:
            delay(1 * ms)
            self.dds.cpld.set_profile(0)  # profile 0 tells CPLD to start ramping up
            delay(2 * ms)
            self.dds.cpld.set_profile(1)  # profile 1 tells CPLD to start ramping back down

