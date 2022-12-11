from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence

from LAX_exp.utilities.conversions import *



class InitializeQubit(LAXSubsequence):
    """
    Subsequence: Initialize Qubit
        Initialize the ion in the S-1/2 mj=-1/2 state.
    """
    name = 'initialize_qubit'

    parameters = {
        'time_spinpol_mu':                  ('timing.time_spinpol_us',                  us_to_mu),
        'time_repump_qubit_mu':             ('timing.time_repump_qubit_us',             us_to_mu),
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu),
        'time_profileswitch_delay_mu':      ('timing.time_profileswitch_delay_us',      us_to_mu)
    }
    devices = [
        'probe',
        'pump',
        'repump_qubit'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # repump pulse
        self.repump_qubit.cfg_sw(1)
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.cfg_sw(0)

        # doppler cooling
        self.pump.cfg_sw(1)
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.cfg_sw(0)

        # spin polarization
        self.probe.cfg_sw(1)
        delay_mu(self.time_spinpol_mu)
        self.probe.cfg_sw(0)
