import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_BIDIR_RAMP, RAM_DEST_ASF, RAM_DEST_POW, RAM_DEST_FTW


class UrukulDitherTest(EnvExperiment):
    """
    Urukul Dither Test.
    Test phase dithering on an AD9910.
    """

    def build(self):
        """
        Ensure that the necessary devices are set.
        """
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch3")

    def prepare(self):
        turns_list_rand = np.random.random(1024)
        self.ram_values = [self.urukul0_ch3.turns_to_pow(rand_val) for rand_val in turns_list_rand]
        self.time_ramp_mu = self.core.seconds_to_mu(2 * us)

    def run(self):
        # program ram
        self._ram_program()

        # set up modulation
        self._ram_setup()

        # start modulating
        self._ram_start()

    @kernel
    def _ram_program(self):
        # set ram profile
        self.urukul0_ch3.set_profile_ram(
            start=0, end=len(self.ram_values)-1, step=1, mode=RAM_MODE_BIDIR_RAMP
        )

        # set profile & update
        self.urukul0_ch3.cpld.set_profile(0)
        self.urukul0_ch3.cpld.io_update.pulse_mu(8)
        delay(1 * ms)

        # write ram
        self.urukul0_ch3.write_ram(self.ram_values)
        delay(100 * ms)

    @kernel
    def _ram_setup(self):
        # write to CFR1 to enable RAM modulation
        self.urukul0_ch3.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_POW)
        self.urukul0_ch3.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set waveform
        self.urukul0_ch3.set_frequency(1 * MHz)
        self.urukul0_ch3.set_att(10 * dB)
        self.urukul0_ch3.cfg_sw(1)
        self.core.break_realtime()

    @kernel
    def _ram_start(self):
        # perform a number of bidirectional ramps
        for i in range(1000):

            # ramp up
            self.core.break_realtime()
            self.urukul0_ch3.cpld.set_profile(0)
            delay_mu(self.time_ramp_mu)

            # ramp down
            self.urukul0_ch3.cpld.set_profile(1)
            self.core.break_realtime()

    def analyze(self):
        pass
