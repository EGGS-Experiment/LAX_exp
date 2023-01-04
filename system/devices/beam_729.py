from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam729(LAXDevice):
    """
    Wrapper for the 729nm qubit beam (polarized).
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "qubit"

    parameters = {
        'freq_qubit_carrier_ftw':           ('beams.freq_mhz.freq_qubit_mhz',               mhz_to_ftw),
        # 'freq_qubit_rsb_ftw':               ('beams.freq_mhz.freq_qubit_rsb_mhz',           mhz_to_ftw),
        # 'freq_qubit_bsb_ftw':               ('beams.freq_mhz.freq_qubit_bsb_mhz',           mhz_to_ftw),
        'ampl_qubit_asf':                   ('beams.ampl_pct.ampl_qubit_pct',               pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul0_ch1'
    }

    @kernel(flags={"fast-math"})
    def prepare_device(self):
        # set carrier profile
        self.core.break_realtime()
        self.beam.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)

        # set RSB profile
        self.core.break_realtime()
        #self.beam.set_mu(self.freq_qubit_rsb_ftw, asf=self.ampl_qubit_asf, profile=1)
        self.beam.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)

        # set BSB profile
        self.core.break_realtime()
        #self.beam.set_mu(self.freq_qubit_bsb_ftw, asf=self.ampl_qubit_asf, profile=2)
        self.beam.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def carrier(self):
        # set carrier profile
        with parallel:
            self.beam.cpld.set_profile(0)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def rsb(self):
        # set red sideband profile
        with parallel:
            self.beam.cpld.set_profile(1)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def bsb(self):
        # set blue sideband profile
        with parallel:
            self.beam.cpld.set_profile(2)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)
