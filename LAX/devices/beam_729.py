from artiq.experiment import *
from LAX_exp.LAX.base_classes import *


class beam_729(Beam_Urukul):
    """
    Wrapper for the 729nm qubit beam (polarized).
    Uses the DDS channel to drive an AOM in double-pass configuration.
    """

    DDS_BOARD =     'beams.dds_board.dds_board_qubit_num'
    DDS_CHANNEL =   'beams.dds_channel.dds_qubit_channel'

    frequencies = [
        "beams.freq_mhz.freq_qubit_carrier_mhz",
        "beams.freq_mhz.freq_qubit_rsb_mhz",
        "beams.freq_mhz.freq_qubit_bsb_mhz",
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_qubit_pct"
    ]


    @kernel(flags='fast-math')
    def _build_set_profiles(self):
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_rsb_ftw, asf=self.ampl_qubit_asf, profile=1)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_bsb_ftw, asf=self.ampl_qubit_asf, profile=2)
