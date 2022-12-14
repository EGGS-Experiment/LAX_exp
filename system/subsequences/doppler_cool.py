from artiq.experiment import *
from LAX_exp.base import LAXSubsequence

from LAX_exp.utilities.conversions import *


class DopplerCool(LAXSubsequence):
    """
    Subsequence: Doppler Cooling
        Cool the ion to the doppler limit using the pump beam.
    """
    name = 'doppler_cool'

    parameters = {
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu)
    }
    devices = [
        'pump'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # doppler cooling
        self.pump.cfg_sw(True)
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.cfg_sw(False)