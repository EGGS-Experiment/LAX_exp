from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


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
        self.pump.on()
        self.pmt.count(self.time_readout_mu)
        self.pump.off()
