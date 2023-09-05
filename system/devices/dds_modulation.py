from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1


class DDSModulation(LAXDevice):
    """
    Device: Modulation DDS

    Uses the DDS channel to modulate the trap RF (used for parametric excitation).
    """
    name = "dds_modulation"
    core_device = ('dds', 'urukul1_ch1')
    devices ={
        'mod_switch': 'ttl11'
    }

    def prepare_device(self):
        self.freq_modulation_ftw =              hz_to_ftw(1 * MHz)
        self.ampl_modulation_asf =              self.get_parameter('ampl_modulation_pct', group='dds.ampl_pct', override=False, conversion_function=pct_to_asf)


    @kernel(flags={"fast-math"})
    def initialize_device(self):
        self.core.break_realtime()

        # close rf switches to kill any modulation signal leakage
        with parallel:
            self.mod_switch.on()
            self.dds.sw.off()

        # set up DDS to reinitialize phase each time we set waveform values
        self.dds.set_mu(self.freq_modulation_ftw, asf=self.ampl_modulation_asf, profile=0)
        self.dds.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.core.break_realtime()

        # todo: set matched latency, cfr1, sine etc.
        # enable matched latency
        self.dds.set_cfr2(matched_latency_enable=1)


    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.dds.sw.on()
            # enable modulation RF switch to DDS
            self.mod_switch.on()

    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.dds.sw.off()
            # disable modulation RF switch for DDS
            self.mod_switch.off()

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.dds.cpld.set_profile(profile_num)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)

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
        # delay_mu()
        # todo: add timing delay for waveform to update
