from artiq.experiment import *
from LAX_exp.LAX.base_classes import *


class beam_854(Beam_Urukul):
    """
    Wrapper for the 854nm qubit repump.
    Uses the DDS channel to drive an AOM.
    """

    DDS_BOARD =     'beams.dds_board.dds_board_num'
    DDS_CHANNEL =   'beams.dds_channel.dds_repump_qubit_channel'

    frequencies = [
        "beams.freq_mhz.freq_repump_qubit_mhz"
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_repump_qubit_pct"
    ]


    @kernel(flags='fast-math')
    def _build_set_profiles(self):
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
