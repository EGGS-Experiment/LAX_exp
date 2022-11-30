from artiq.experiment import *
from LAX_exp.LAX.base_classes import *


class beam_397_probe(Beam_Urukul):
    """
    Wrapper for the 397nm probe beam (polarized).
    Uses the DDS channel to drive an AOM in double-pass configuration.
    """

    DDS_BOARD =     'beams.dds_board.dds_board_num'
    DDS_CHANNEL =   'beams.dds_channel.dds_probe_channel'

    frequencies = [
        "beams.freq_mhz.freq_probe_redist_mhz"
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_probe_redist_pct"
    ]


    @kernel(flags='fast-math')
    def _build_set_profiles(self):
        self.core.break_realtime()
        self.dev.set_mu(self.freq_probe_redist_ftw, asf=self.ampl_probe_redist_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_probe_redist_ftw, asf=self.ampl_probe_redist_asf, profile=1)
