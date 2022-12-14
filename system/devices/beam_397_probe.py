from artiq.experiment import *

from LAX_exp.base import LAXDevice
from LAX_exp.utilities.conversions import *


class Beam397Probe(LAXDevice):
    """
    Wrapper for the 397nm probe beam (polarized).
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "probe"

    parameters = {
        'freq_spinpol_ftw':          ('beams.freq_mhz.freq_probe_spinpol_mhz',      mhz_to_ftw),
        'ampl_spinpol_asf':          ('beams.ampl_pct.ampl_probe_spinpol_pct',      pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch0'
    }

    @kernel(flags='fast-math')
    def prepare_device(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=1)

    @kernel(flags='fast-math')
    def on(self):
        self.beam.cfg_sw(1)

    @kernel(flags='fast-math')
    def off(self):
        self.beam.cfg_sw(0)
