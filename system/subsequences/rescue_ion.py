from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class RescueIon(LAXSubsequence):
    """
    Subsequence: Rescue Ion

    Rescue the ion by running the 397nm cooling beam at rescue parameters.
    """
    name = 'rescue_ion'

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')

    def prepare_subsequence(self):
        self.time_rescue_mu = self.get_parameter('time_rescue_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self):
        # set rescue waveform and ensure 866 is on
        self.pump.rescue()
        self.repump_cooling.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_rescue_mu)
        self.pump.off()
