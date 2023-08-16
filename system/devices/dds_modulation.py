from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE


class DDSModulation(LAXDevice):
    """
    Device: Modulation DDS

    Uses the DDS channel to modulate the trap RF (used for parametric excitation).
    """
    name = "dds_modulation"
    core_device = ('dds', 'urukul0_ch2')
    devices ={
        'mod_switch': 'ttl12'
    }

    def prepare_device(self):
        self.freq_modulation_ftw =              hz_to_ftw(1 * MHz)
        self.ampl_modulation_asf =              self.get_parameter('ampl_modulation_pct', group='dds.ampl_pct', override=False, conversion_function=pct_to_asf)


    @kernel(flags={"fast-math"})
    def initialize_device(self):
        self.core.break_realtime()

        # enable modulation signal to feed through
        with parallel:
            self.mod_switch.on()
            self.cfg_sw(False)

        # set up DDS to reinitialize phase each time we set waveform values
        self.set_mu(self.freq_modulation_ftw, asf=self.ampl_modulation_asf, profile=0)
        self.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.cfg_sw(True)
            # enable modulation RF switch to DDS
            self.mod_switch.on()

    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.cfg_sw(False)
            # disable modulation RF switch for DDS
            self.mod_switch.off()

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.dds.cpld.set_profile(profile_num)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def reset_phase(self):
        """
        todo: document
        :return:
        """
        self.set_cfr1(phase_autoclear=1)
        self.dds.cpld.io_update.pulse_mu(8)
