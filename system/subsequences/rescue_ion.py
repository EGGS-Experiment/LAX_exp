from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class RescueIon(LAXSubsequence):
    """
    Subsequence: Rescue Ion

    Rescue the ion by running the 397nm cooling beam at rescue parameters.
    """
    name = 'rescue_ion'

    parameters = {
        'time_rescue_mu':                   ('timing.time_rescue_us',          us_to_mu)
    }

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('cooling_repump')

    @kernel(flags={"fast-math"})
    def run(self):
        # set rescue waveform and ensure 866 is on
        self.pump.rescue()
        self.repump.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_rescue_mu)
        self.pump.off()
