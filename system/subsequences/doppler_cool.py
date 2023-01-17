from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class DopplerCool(LAXSubsequence):
    """
    Subsequence: Doppler Cooling
        Cool the ion to the doppler limit using the pump beam.
    """
    name = 'doppler_cool'

    parameters = {
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu)
    }

    def build_subsequence(self):
        self.setattr_device('pump')

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()
