from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class InitializeQubit(LAXSubsequence):
    """
    Subsequence: Initialize Qubit
        Initialize the ion in the S-1/2 mj=-1/2 state.
    """
    name = 'initialize_qubit'

    parameters = {
        'time_spinpol_mu':                  ('timing.time_spinpol_us',                  us_to_mu),
        'time_repump_qubit_mu':             ('timing.time_repump_qubit_us',             us_to_mu),
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu)
    }
    devices = [
        'probe',
        'pump',
        'repump_cooling',
        'repump_qubit'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # permanently enable cooling repump
        self.repump_cooling.on()

        # repump pulse
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.off()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # spin polarization
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
