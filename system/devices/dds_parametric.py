from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class DDSParametric(LAXDevice):
    """
    Device: Parametric DDS

    Uses the DDS channel to modulate the trap RF (used for parametric excitation).
    """
    name = "dds_parametric"
    core_device = ('dds', 'urukul1_ch1')
    devices = {
        'mod_switch': 'ttl11'
    }
    kernel_invariants = {
        "cpld", "sw",
        "ampl_modulation_asf", "freq_cleanup_ftw"
    }

    def prepare_device(self):
        # break out AD9910 attributes/devices
        self.sw =   self.dds.sw
        self.cpld = self.dds.cpld

        # get relevant parameters
        self.ampl_modulation_asf =  self.get_parameter('ampl_parametric_pct', group='dds.ampl_pct',
                                                       override=False, conversion_function=pct_to_asf)
        self.freq_cleanup_ftw =     self.dds.frequency_to_ftw(150 * MHz)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # close rf switches to kill any modulation signal leakage
        self.mod_switch.off()
        self.sw.off()

        # ensure output has matched latency
        self.dds.set_cfr2(matched_latency_enable=1)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        # set default profile
        self.set_profile(DEFAULT_PROFILE)
        delay_mu(8000)

        # clear any possible output
        self.dds.set_att_mu(0)
        self.dds.set_mu(self.freq_cleanup_ftw, asf=0x01, profile=DEFAULT_PROFILE)
        delay_mu(10000)

        # make sure switches are closed
        self.off()
        delay_mu(5000)

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        with parallel:
            # enable RF switch onboard Urukul
            self.sw.on()
            # enable modulation RF switch to DDS
            self.mod_switch.on()

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        with parallel:
            # disable RF switch onboard Urukul
            self.sw.off()
            # disable modulation RF switch for DDS
            self.mod_switch.off()

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        self.cpld.set_profile(profile_num)
        self.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_phase_absolute(self) -> TNone:
        """
        todo: document
        """
        # set up DDS to reinitialize phase each time we set waveform values
        self.dds.set_phase_mode(PHASE_MODE_ABSOLUTE)

        # enable matched latency and phase autoclearing
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13) |    # phase_autoclear
                         2)             # default serial I/O configs

    @kernel(flags={"fast-math"})
    def reset_phase(self) -> TNone:
        """
        todo: document
        """
        # ensure signal is output as a sine with 0 phase
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13) |    # phase_autoclear
                         2)             # default serial I/O configs
        # pulse io_update to clear phase
        at_mu(now_mu() & ~0x7)
        self.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_AD9910_PHASE_AUTOCLEAR_DELAY_MU)

    @kernel(flags={"fast-math"})
    def io_update(self) -> TNone:
        """
        Pulse the CPLDs IO_UPDATE pin.
        Can be used to clear the phase accumulator if the phase_autoclear
            flag is set in CFR1.
        """
        self.cpld.io_update.pulse_mu(8)
