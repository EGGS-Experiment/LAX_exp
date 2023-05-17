"""
Contains specifically constructed subsequences for LAX.
"""

__all__ = []

# general
from LAX_exp.system.subsequences.rescue_ion import RescueIon
from LAX_exp.system.subsequences.cleanup import Cleanup
from LAX_exp.system.subsequences.no_operation import NoOperation
__all__.extend(['RescueIon', 'Cleanup', 'NoOperation'])

# motional
from LAX_exp.system.subsequences.parametric_excite import ParametricExcite
__all__.extend(['ParametricExcite'])

# cooling
from LAX_exp.system.subsequences.doppler_cool import DopplerCool
from LAX_exp.system.subsequences.sideband_cool_pulsed import SidebandCoolPulsed
from LAX_exp.system.subsequences.sideband_cool_continuous import SidebandCoolContinuous
__all__.extend(['DopplerCool', 'SidebandCoolPulsed', 'SidebandCoolContinuous'])

# state preparation
from LAX_exp.system.subsequences.spin_polarization import SpinPolarization
from LAX_exp.system.subsequences.initialize_qubit import InitializeQubit
__all__.extend(['SpinPolarization', 'InitializeQubit'])

# state manipulation
from LAX_exp.system.subsequences.tickle import Tickle
from LAX_exp.system.subsequences.rabi_flop import RabiFlop
from LAX_exp.system.subsequences.ramsey import Ramsey
# from LAX_exp.system.subsequences.phaser_pulse import PhaserPulse
__all__.extend(['Tickle', 'RabiFlop', 'Ramsey'])

# readout
from LAX_exp.system.subsequences.readout import Readout
from LAX_exp.system.subsequences.absorption_probe import AbsorptionProbe
from LAX_exp.system.subsequences.absorption_probe2 import AbsorptionProbe2
__all__.extend(['AbsorptionProbe', 'AbsorptionProbe2'])
