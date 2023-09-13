from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class BeamManager(LAXDevice):
    """
    Device: Beam Manager

    A meta-device that manages the beam DDSs relevant to the ion.
    Used to reduce overhead and code complexity, making it more straightforward to
    run composite pulse sequences.
    """
    name = "qubit"
    core_device = ('beam', 'urukul2_ch1')
    devices ={
        'rf_switch':    'ttl12'
    }

    def prepare_device(self):
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.beam.cfg_sw(True)

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.beam.cfg_sw(False)

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.beam.cpld.set_profile(profile_num)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
