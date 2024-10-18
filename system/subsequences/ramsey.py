from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Ramsey(LAXSubsequence):
    """
    Subsequence: Ramsey

    Do two pi/2 pulses, separated by a given delay.
    """
    name = 'ramsey'
    kernel_invariants = {
        "time_pi2pulse_mu", "time_ramsey_delay_mu"
    }

    def build_subsequence(self):
        # get devices
        self.setattr_device('qubit')

        # get arguments
        self.setattr_argument('time_pi2pulse_us',       NumberValue(default=125, precision=3, step=10, min=1, max=1000000), group='ramsey_spectroscopy')
        self.setattr_argument('time_ramsey_delay_us',   NumberValue(default=1000, precision=3, step=10, min=1, max=1000000), group='ramsey_spectroscopy')

    def prepare_subsequence(self):
        # convert parameters to machine units
        self.time_pi2pulse_mu =     self.core.seconds_to_mu(self.time_pi2pulse_us * us)
        self.time_ramsey_delay_mu = self.core.seconds_to_mu(self.time_ramsey_delay_us * us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # initial pi/2 pulse
        self.qubit.on()
        delay_mu(self.time_pi2pulse_us)
        self.qubit.off()

        # ramsey delay
        delay_mu(self.time_ramsey_delay_mu)

        # final pi/2 pulse
        self.qubit.on()
        delay_mu(self.time_half_pipulse_mu)
        self.qubit.off()
