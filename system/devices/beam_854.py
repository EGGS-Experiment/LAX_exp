from artiq.experiment import *

from LAX_exp.utilities.conversions import *
from LAX_exp.base import LAXDevice


class Beam854(LAXDevice):
    """
    Wrapper for the 854nm qubit repump.
        Uses the DDS channel to drive an AOM.
    """
    name = "qubit_repump"

    parameters = {
        'freq_repump_qubit_ftw':            ('beams.freq_mhz.freq_repump_qubit_mhz',        mhz_to_ftw),
        'ampl_repump_qubit_asf':            ('beams.ampl_pct.ampl_repump_qubit_pct',        pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch3'
    }

    @kernel(flags='fast-math')
    def prepare_device(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)

    @kernel(flags='fast-math')
    def on(self):
        self.beam.cfg_sw(1)

    @kernel(flags='fast-math')
    def off(self):
        self.beam.cfg_sw(0)
