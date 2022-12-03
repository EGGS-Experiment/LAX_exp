from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class Readout(LAXSubsequence):
    """
    Subsequence: Readout
        Apply the 397nm pump beam while reading fluorescence via the PMT.
    """
    name = 'readout'

    devices = [
        'pump',
        'pmt'
    ]
    subsequence_parameters = {
        'time_readout_mu':              ('timing.time_readout_us', us_to_mu),
        'time_profileswitch_delay_mu':  ('timing.time_profileswitch_delay_us', us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform
        with parallel:
            self.dds_board.set_profile(1)
            delay_mu(self.time_profileswitch_delay_mu)

        # readout pulse
        self.pump.cfg_sw(1)
        self.pmt.count(self.time_readout_mu)
        self.pump.cfg_sw(0)
