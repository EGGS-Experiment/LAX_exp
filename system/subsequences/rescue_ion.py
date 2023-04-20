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
        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')

        # rescue cooling (happens at end of each repetition)
        self.setattr_argument("repetitions_per_cooling",                NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000), group='rescue_ion')
        self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000), group='rescue_ion')

    def prepare_subsequence(self):
        self.time_rescue_mu =                                           self.get_parameter('time_rescue_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self, i):
        # check whether it's time to rescue the ion
        if (i > 0) and (i % self.repetitions_per_cooling == 0):

            # set rescue waveform and ensure 866 is on
            self.pump.rescue()
            self.repump_cooling.on()

            # doppler cooling
            self.pump.on()
            delay_mu(self.time_rescue_mu)
            self.pump.off()
