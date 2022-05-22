import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_BIDIR_RAMP, RAM_DEST_ASF


class Urukul_AWG(EnvExperiment):
    """
    Urukul AWG Test.
    """

    def build(self):
        """
        Ensure that the necessary devices are set.
        """
        self.setattr_device("core")                                 # always needed
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")

    def prepare(self):
        # todo: get amplitude RAM data
        self.ram_values = np.linspace(0, 0x3FFF, 0x2FF, dtype=np.int32)

    @kernel
    def run(self):
        self.core.reset()

        # initialize devices
        self.urukul0_cpld.init()
        self.core.break_realtime()
        self.urukul0_ch0.init()
        self.core.break_realtime()

        # set ram profile
        self.urukul0_ch0.set_profile_ram(
            start=0, end=len(self.ram_values)-1, step=0xFFFFFFFF, mode=RAM_MODE_BIDIR_RAMP
        )

        # set profile & update
        self.urukul0_ch0.cpld.set_profile(0)
        self.urukul0_ch0.io_update.pulse(20 * ns)
        self.core.break_realtime()

        # write ram
        self.urukul0_ch0.write_ram(self.ram_values)
        delay(10 * ms)

        # write to CFR1 to enable RAM modulation
        self.urukul0_ch0.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
        self.urukul0_ch0.cpld.io_update.pulse(20 * ns)
        self.core.break_realtime()

        # set waveform
        self.set_frequency(100 * MHz)
        self.set_att(10 * dB)
        self.urukul0_ch0.cfg_sw(1)
        self.core.break_realtime()

    def analyze(self):
        pass
