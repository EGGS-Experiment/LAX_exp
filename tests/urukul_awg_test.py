import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_BIDIR_RAMP, RAM_DEST_ASF


class Urukul_AWG(EnvExperiment):
    """
    Urukul AWG Test
    Test the use of the Urukul as an AWG.
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
        #self.ram_values = np.linspace(0x0FF, 0x3FFF, 1024, dtype=np.int32)
        pass

    @kernel
    def run(self):
        self.core.reset()

        # create data
        ram_values = [0] * 1024
        for i in range(len(ram_values)):
            ram_values[i] = 13 * i

        self.core.break_realtime()

        # initialize devices
        self.urukul0_cpld.init()
        self.core.break_realtime()
        self.urukul0_ch0.init()
        self.core.break_realtime()

        # set ram profile
        self.urukul0_ch0.set_profile_ram(
            start=0, end=len(ram_values)-1, step=1, mode=RAM_MODE_BIDIR_RAMP
        )

        # set profile & update
        self.urukul0_ch0.cpld.set_profile(0)
        self.urukul0_ch0.cpld.io_update.pulse_mu(8)
        delay(1 * ms)

        # write ram
        self.urukul0_ch0.write_ram(ram_values)
        delay(100 * ms)

        # write to CFR1 to enable RAM modulation
        self.urukul0_ch0.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
        self.urukul0_ch0.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set waveform
        self.urukul0_ch0.set_frequency(100 * MHz)
        self.urukul0_ch0.set_att(5 * dB)
        self.urukul0_ch0.cfg_sw(1)
        self.core.break_realtime()

        while True:
            self.core.break_realtime()
            self.urukul0_ch0.cpld.set_profile(0)
            delay(2 * us)
            self.urukul0_ch0.cpld.set_profile(1)
            self.core.break_realtime()

    def analyze(self):
        pass
