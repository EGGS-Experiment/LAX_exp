from artiq.experiment import *
from LAX_exp.LAX.base_classes import *
# todo: maybe set profile names as enum (e.g. set_profile(PROFILE_COOLING))
# gate counts
# get counts
# read out
# todo: core_dma must only happen at sequence level and above

class beam_397_pump(Beam_Urukul):
    """
    A beam
    # todo: document
    """

    DDS_NAME = 'urukul1_ch1'

    frequencies = [
        "beams.freq_mhz.freq_pump_cooling_mhz",
        "beams.freq_mhz.freq_pump_readout_mhz"
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_pump_cooling_pct",
        "beams.ampl_pct.ampl_pump_readout_pct"
    ]

    @kernel
    def _build_set_profiles(self):
        self.dev.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dev.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()
