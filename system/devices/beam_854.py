from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam854(LAXDevice):
    """
    Device: Qubit Repump Beam (854nm)

    Uses the DDS channel to drive an AOM.
    """
    name = "repump_qubit"
    core_device = ('beam', 'urukul2_ch3')
    devices = {
        'rf_switch':    'ttl13'
    }
    kernel_invariants = {
        "freq_repump_qubit_ftw",
        "ampl_repump_qubit_asf"
    }

    def prepare_device(self):
        self.freq_repump_qubit_ftw = self.get_parameter('freq_repump_qubit_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_repump_qubit_asf = self.get_parameter('ampl_repump_qubit_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=2)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        with parallel:
            # enable RF switch onboard Urukul
            self.beam.sw.on()

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        with parallel:
            # disable RF switch onboard Urukul
            self.beam.sw.off()

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        self.beam.cpld.set_profile(profile_num)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
