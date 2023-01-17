from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam866(LAXDevice):
    """
    Wrapper for the 866nm cooling repump.
        Uses the DDS channel to drive an AOM.
    """
    name = "cooling_repump"

    parameters = {
        'freq_repump_cooling_ftw':          ('beams.freq_mhz.freq_repump_cooling_mhz',      mhz_to_ftw),
        'ampl_repump_cooling_asf':          ('beams.ampl_pct.ampl_repump_cooling_pct',      pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch2'
    }

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)
