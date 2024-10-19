from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Readout(LAXSubsequence):
    """
    Subsequence: Readout

    Read out the ion state by shining the pump beam and reading fluorescence via PMT counts.
    """
    name = 'readout'
    kernel_invariants = {
        "time_readout_mu"
    }

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        self.time_readout_mu = self.get_parameter('time_readout_us', group='timing',
                                                  override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # set readout waveform
        self.pump.readout()

        # readout pulse
        self.pump.on()
        self.repump_cooling.on()
        self.pmt.count(self.time_readout_mu)
        # tmp remove - leave 397nm readout on for max cooling
        # self.pump.off()

    @kernel(flags={"fast-math"})
    def fetch_count(self) -> TInt32:
        """
        Convenience function so that users don't have to separately instantiate the PMT
        device object to read counts.
        """
        return self.pmt.fetch_count()
