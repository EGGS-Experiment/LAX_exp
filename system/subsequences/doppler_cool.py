from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class DopplerCool(LAXSubsequence):
    """
    Subsequence: Doppler Cooling

    Cool the ion to the doppler limit using the pump beam.
    """
    name = 'doppler_cool'
    kernel_invariants = {
        "time_doppler_cooling_mu"
    }

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

    def prepare_subsequence(self):
        # get doppler cooling time
        self.time_doppler_cooling_mu = self.get_parameter('time_doppler_cooling_us', group='timing',
                                                          override=False, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # turn repumps on
        self.repump_cooling.on()
        self.repump_qubit.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()
