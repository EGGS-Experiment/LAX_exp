from artiq.experiment import *
from inspect import getmembers, ismethod

from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf, us_to_mu


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
        'ampl_qubit_asf':                   ('beams.ampl_pct.ampl_qubit_pct',               pct_to_asf),
        'time_profileswitch_delay_mu':      ('timing.time_profileswitch_delay_us',          us_to_mu)

    }
    core_devices = {
        'beam': 'urukul0_ch2'
    }

    @kernel(flags='fast-math')
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

    @kernel(flags='fast-math')
    def on(self):
        self.beam.cfg_sw(1)

    @kernel(flags='fast-math')
    def off(self):
        self.beam.cfg_sw(0)

    @kernel(flags='fast-math')
    def carrier(self):
        # set carrier profile
        with parallel:
            self.beam.cpld.set_profile(0)
            delay_mu(self.time_profileswitch_delay_mu)

    @kernel(flags='fast-math')
    def rsb(self):
        # set red sideband profile
        with parallel:
            self.beam.cpld.set_profile(1)
            delay_mu(self.time_profileswitch_delay_mu)

    @kernel(flags='fast-math')
    def bsb(self):
        # set blue sideband profile
        with parallel:
            self.beam.cpld.set_profile(2)
            delay_mu(self.time_profileswitch_delay_mu)
