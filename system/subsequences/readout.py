from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Readout(LAXSubsequence):
    """
    Subsequence: Readout
        Read out the ion state by shining the pump beam and reading fluorescence via PMT counts.
    """
    name = 'readout'

    devices = [
        'pump',
        'pmt'
    ]
    parameters = {
        'time_readout_mu':                  ('timing.time_readout_us',                  us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform
        self.pump.readout()

        # readout pulse
        self.pump.cfg_sw(True)
        self.pmt.count(self.time_readout_mu)
        self.pump.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def fetch_count(self):
        # convenience function so that users don't have to separately instantiate the PMT
        # device object to read counts
        return self.pmt.fetch_count()
