from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam729(LAXDevice):
    """
    Device: Qubit Beam (729nm, polarized)

    Uses the DDS channel to drive the 729nm AOM in double-pass configuration.
    """
    name = "qubit"
    core_device = ('beam', 'urukul0_ch0')
    devices = {
        'rf_switch':    'ttl14'
    }
    kernel_invariants = {
        "freq_qubit_ftw",
        "ampl_qubit_asf"
    }

    def prepare_device(self):
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.beam.sw.on()

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.beam.sw.off()

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.beam.cpld.set_profile(profile_num)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
