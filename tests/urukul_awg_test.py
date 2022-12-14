import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_BIDIR_RAMP, RAM_DEST_ASF, RAM_DEST_FTW


class Urukul_AWG(EnvExperiment):
    """
    Urukul AWG Test
    Test the use of the Urukul as an AWG.
    """

    def build(self):
        """
        Ensure that the necessary devices are set.
        """
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.dds = self.get_device("urukul0_ch2")

    def prepare(self):
        # get ram data
        #self.ram_values = np.linspace(0x0FF, 0x3FFF, 1024, dtype=np.int32)
        self.ram_values = np.linspace(0x1C28F5C2, 0x1D70A3D6, 1024, dtype=np.int32)
        pass

    @kernel
    def run(self):
        self.core.reset()

        # set ram profile
        self.dds.set_profile_ram(
            start=0,
            end=len(self.ram_values)-1,
            step=1,
            mode=RAM_MODE_BIDIR_RAMP
        )

        # set profile & update
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay(10 * ms)

        # write ram
        self.dds.write_ram(self.ram_values)
        delay(100 * ms)

        # write to CFR1 to enable RAM modulation
        self.dds.set_cfr1(
            ram_enable=1,
            ram_destination=RAM_DEST_FTW
        )
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set waveform
        self.dds.set_frequency(110 * MHz)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_att(8 * dB)
        self.dds.cfg_sw(True)
        self.core.break_realtime()

        while True:
            delay(1 * ms)
            self.dds.cpld.set_profile(0)
            delay(1 * ms)
            self.dds.cpld.set_profile(1)
            self.core.break_realtime()
