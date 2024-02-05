from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1


class DDSDipole(LAXDevice):
    """
    Device: Dipole Modulation DDS

    Uses the DDS channel to apply a dipole field on the EGGS-side.
    """
    name = "dds_dipole"
    core_device = ('dds', 'urukul1_ch3')
    kernel_invariants = {
        "freq_modulation_ftw",
        "ampl_modulation_asf"
    }

    def prepare_device(self):
        self.freq_modulation_ftw =              hz_to_ftw(1 * MHz)
        self.ampl_modulation_asf =              self.get_parameter('ampl_modulation_pct', group='dds.ampl_pct',
                                                                   override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        self.core.break_realtime()

        # close rf switches to kill any modulation signal leakage
        self.dds.sw.off()

        # set up DDS to reinitialize phase each time we set waveform values
        self.dds.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.dds.set_mu(self.freq_modulation_ftw, asf=self.ampl_modulation_asf, profile=0)
        self.core.break_realtime()

        # todo: set matched latency, cfr1, sine etc.
        # enable matched latency
        self.dds.set_cfr2(matched_latency_enable=1)


    @kernel(flags={"fast-math"})
    def on(self):
        self.dds.sw.on()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def off(self):
        self.dds.sw.off()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.dds.cpld.set_profile(profile_num)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def reset_phase(self):
        """
        todo: document
        :return:
        """
        # ensure signal is output as a sine with 0 phase
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13))     # phase_autoclear
        self.dds.cpld.io_update.pulse_mu(8)
