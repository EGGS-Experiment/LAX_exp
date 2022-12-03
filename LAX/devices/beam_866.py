from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam866(LAXDevice):
    """
    Wrapper for the 866nm cooling repump.
        Uses the DDS channel to drive an AOM.
    """

    device_names = {'beam': 'urukul1_ch2'}
    device_parameters = {
        'freq_repump_cooling_ftw': ('beams.freq_mhz.freq_repump_cooling_mhz', mhz_to_ftw),
        'ampl_repump_cooling_saf': ('beams.ampl_pct.ampl_repump_cooling_pct', pct_to_asf)
    }

    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
