from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam854(LAXDevice):
    """
    Wrapper for the 854nm qubit repump.
        Uses the DDS channel to drive an AOM.
    """
    name = "qubit_repump"

    device_names = {'beam': 'urukul1_ch3'}
    device_parameters = {
        'freq_repump_qubit_ftw': ('beams.freq_mhz.freq_repump_qubit_mhz', mhz_to_ftw),
        'ampl_repump_qubit_asf': ('beams.ampl_pct.ampl_repump_qubit_pct', pct_to_asf)
    }


    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
