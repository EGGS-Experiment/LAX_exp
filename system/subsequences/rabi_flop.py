from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class RabiFlop(LAXSubsequence):
    """
    Subsequence: Rabi Flop

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    """
    name = 'rabi_flop'

    parameters = {
        'time_rabiflop_mu':                 ('timing.time_rabiflop_us',                 us_to_mu)
    }
    devices = [
        'qubit'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform
        # todo: remove this line
        #self.qubit.carrier()

        # population transfer pulse
        self.qubit.on()
        #delay_mu(self.time_rabiflop_mu)
        delay_mu(400000)
        self.qubit.off()
