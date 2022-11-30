from artiq.experiment import *
from LAX_exp.LAX.base_classes import *


class beam_866(Beam_Urukul):
    """
    Wrapper for the 866nm cooling repump.
    Uses the DDS channel to drive an AOM.
    """

    DDS_BOARD =     'beams.dds_board.dds_board_num'
    DDS_CHANNEL =   'beams.dds_channel.dds_repump_cooling_channel'

    frequencies = [
        "beams.freq_mhz.freq_repump_cooling_mhz"
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_repump_cooling_pct"
    ]


    @kernel(flags='fast-math')
    def _build_set_profiles(self):
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
