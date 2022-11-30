from artiq.experiment import *
from LAX_exp.LAX.base_classes import *
# todo: maybe set profile names as enum (e.g. set_profile(PROFILE_COOLING))


class beam_397_pump(Beam_Urukul):
    """
    Wrapper for the 397nm pump beam.
    Uses the DDS channel to drive an AOM in double-pass configuration.
    """

    DDS_BOARD =     'beams.dds_board.dds_board_num'
    DDS_CHANNEL =   'beams.dds_channel.dds_pump_channel'

    frequencies = [
        "beams.freq_mhz.freq_pump_cooling_mhz",
        "beams.freq_mhz.freq_pump_readout_mhz"
    ]

    amplitudes = [
        "beams.ampl_pct.ampl_pump_cooling_pct",
        "beams.ampl_pct.ampl_pump_readout_pct"
    ]


    @kernel(flags='fast-math')
    def _build_set_profiles(self):
        self.core.break_realtime()
        self.dev.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
