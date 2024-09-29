from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1


class DDSParametric(LAXDevice):
    """
    Device: Parametric DDS

    Uses the DDS channel to modulate the trap RF (used for parametric excitation).
    """
    name = "dds_parametric"
    core_device = ('dds', 'urukul1_ch1')
    devices = {
        'mod_switch': 'ttl11',
        'servo_hold': 'ttl10'
    }
    kernel_invariants = {
        "ampl_modulation_asf",
        "sw"
    }

    def prepare_device(self):
        self.ampl_modulation_asf = self.get_parameter('ampl_parametric_pct', group='dds.ampl_pct', override=False, conversion_function=pct_to_asf)
        # break out AD9910 attributes/devices
        self.sw = self.dds.sw
        self.cpld = self.dds.cpld

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # close rf switches to kill any modulation signal leakage
        self.mod_switch.off()
        self.dds.sw.off()

        # ensure output has matched latency
        self.dds.set_cfr2(matched_latency_enable=1)


    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.dds.sw.on()
            # enable modulation RF switch to DDS
            self.mod_switch.on()
            # enable integrator hold for the trap RF servo
            # self.servo_hold.on()

    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.dds.sw.off()
            # disable modulation RF switch for DDS
            self.mod_switch.off()
            # resume trap RF servo
            # self.servo_hold.off()

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32):
        self.dds.cpld.set_profile(profile_num)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_phase_absolute(self):
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
    def reset_phase(self):
        """
        todo: document
        :return:
        """
        # ensure signal is output as a sine with 0 phase
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 16) |    # select_sine_output
                         (1 << 13))     # phase_autoclear
        # pulse io_update to clear phase
        at_mu(now_mu() & ~0x7)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_AD9910_PHASE_AUTOCLEAR_DELAY_MU)

    @kernel(flags={"fast-math"})
    def io_update(self):
        """
        Pulse the CPLDs IO_UPDATE pin.
        Can be used to clear the phase accumulator if the phase_autoclear
            flag is set in CFR1.
        """
        self.dds.cpld.io_update.pulse_mu(8)
