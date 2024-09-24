from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1


class Ramsey(LAXSubsequence):
    """
    Subsequence: Ramsey

    Do two pi/2 pulses, separated by a given delay.
    """
    name = 'ramsey'
    kernel_invariants = {
        "time_ramsey_pulse_mu",
        "time_ramsey_delay_mu",
        "phase_ramsey_pow"
    }

    def build_subsequence(self):
        # get arguments
        self.setattr_argument('time_ramsey_pulse_us',       NumberValue(default=5, ndecimals=3, step=10, min=4, max=1000000), group='ramsey_spectroscopy')
        self.setattr_argument('time_ramsey_delay_us',       NumberValue(default=125, ndecimals=3, step=10, min=1, max=1000000), group='ramsey_spectroscopy')
        self.setattr_argument('phase_ramsey_turns',         NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ramsey_spectroscopy')

        # get devices
        self.setattr_device('qubit')

    def prepare_subsequence(self):
        # convert parameters to machine units
        self.time_ramsey_pulse_mu =     self.core.seconds_to_mu(self.time_ramsey_pulse_us * us)
        self.time_ramsey_delay_mu =     self.core.seconds_to_mu(self.time_ramsey_delay_us * us)
        self.phase_ramsey_pow =         self.qubit.turns_to_pow(self.phase_ramsey_turns)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # ensure qubit beam has matched latencies
        self.qubit.set_cfr2(matched_latency_enable=1)

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the ramsey pulse sequence.
        Ensures phase coherence between the pulses.
        """
        '''
        PREPARE WAVEFORM
        '''
        # ensure phase_autoclear is enabled ahead of time
        self.qubit.write32(_AD9910_REG_CFR1,
                           (1 << 16) |  # select_sine_output
                           (1 << 13))   # phase_autoclear

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7
        # begin output waveform
        at_mu(time_start_mu)
        self.qubit.set_profile(0)

         #todo: delay sufficient time


        '''
        RAMSEY SEQUENCE
        '''
        # initial pulse
        self.qubit.on()
        delay_mu(self.time_ramsey_pulse_mu)
        self.qubit.off()

        # ramsey delay
        delay_mu(self.time_ramsey_delay_mu)

         # todo: add phase delay

        # reverse pulse
        self.qubit.on()
        delay_mu(self.time_ramsey_pulse_mu)
        self.qubit.off()
