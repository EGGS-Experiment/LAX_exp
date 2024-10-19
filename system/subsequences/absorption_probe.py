from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

# tmp remove
import numpy as np
# tmp remove


class AbsorptionProbe(LAXSubsequence):
    """
    Subsequence: Absorption Probe

    Cool the ion before switching to a short, low-power pulse to measure absorption on a transition.
    """
    name = 'absorption_probe'
    kernel_invariants = {
        "time_doppler_cooling_mu",
        "time_probe_mu"
    }

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        self.time_doppler_cooling_mu =      self.get_parameter('time_doppler_cooling_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_probe_mu =                self.get_parameter('time_probe_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # set cooling waveform
        self.pump.cooling()
        # enable 854nm during sequence to evacuate D-5/2 population
        self.repump_qubit.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # switch to probe waveform
        self.pump.set_profile(1)

        # record fluorescence
        self.pump.on()
        self.pmt.count(self.time_probe_mu)
        self.pump.off()
