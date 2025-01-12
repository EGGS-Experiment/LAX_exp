from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class DDSDipole(LAXDevice):
    """
    Device: Dipole Modulation DDS

    Uses the DDS channel to apply a dipole field on the EGGS-side.
    """
    name = "dds_dipole"
    core_device = ('dds', 'urukul1_ch3')
    devices = {
        'rf_switch': 'ttl15'
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
        self.ampl_modulation_asf =  self.get_parameter('ampl_dipole_pct', group='dds.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.freq_cleanup_ftw =     self.dds.frequency_to_ftw(100 * MHz)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # close rf switches to kill any modulation signal leakage
        self.sw.off()

        # set up DDS to reinitialize phase each time we set waveform values
        self.dds.set_phase_mode(PHASE_MODE_ABSOLUTE)

        # enable matched latency and phase autoclearing
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13))     # phase_autoclear
        self.dds.set_cfr2(matched_latency_enable=1)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        self.core.break_realtime()

        # set default profile
        self.set_profile(DEFAULT_PROFILE)
        self.cpld.io_update.pulse_mu(8)

        # clear any possible output
        self.dds.set_att_mu(0x0)
        self.core.break_realtime()
        self.dds.set_mu(self.freq_cleanup_ftw, asf=0x01, profile=DEFAULT_PROFILE)
        self.core.break_realtime()

        # make sure switches are closed
        self.off()

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        # note: didn't bother parallel-ing this section since it's unimportant
        # enable RF switch onboard Urukul
        self.sw.on()
        delay_mu(1)

        # enable external RF switch
        self.rf_switch.on()
        delay_mu(20000)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        # note: didn't bother parallel-ing this section since it's unimportant
        # enable RF switch onboard Urukul
        self.sw.off()
        delay_mu(1)

        # enable external RF switch
        self.rf_switch.off()
        delay_mu(20000)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        self.cpld.set_profile(profile_num)
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
                         (1 << 13))     # phase_autoclear

    @kernel(flags={"fast-math"})
    def reset_phase(self) -> TNone:
        """
        todo: document
        """
        # ensure signal is output as a sine with 0 phase
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13))     # phase_autoclear

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
