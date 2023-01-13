from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam729(LAXDevice):
    """
    Device: Wrapper for the 729nm qubit beam (polarized).

    Uses the DDS channel to drive the 729nm AOM in double-pass configuration.
    """
    name = "qubit"

    parameters = {
        'freq_qubit_carrier_ftw':           ('beams.freq_mhz.freq_qubit_mhz',               mhz_to_ftw),
        'ampl_qubit_asf':                   ('beams.ampl_pct.ampl_qubit_pct',               pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul0_ch1'
    }

    @kernel(flags={"fast-math"})
    def prepare_device(self):
        self.core.break_realtime()
        self.beam.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def carrier(self):
        self.core.break_realtime()

        # set carrier profile
        with parallel:
            self.beam.cpld.set_profile(0)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)

        self.core.break_realtime()
