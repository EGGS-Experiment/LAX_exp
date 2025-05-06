"""
Contains specifically constructed subsequences for LAX.
"""

__all__ = []

# general
from LAX_exp.system.subsequences.rescue_ion import RescueIon
from LAX_exp.system.subsequences.cleanup import Cleanup
from LAX_exp.system.subsequences.no_operation import NoOperation
__all__.extend(['RescueIon', 'Cleanup', 'NoOperation'])

# cooling
from LAX_exp.system.subsequences.doppler_cool import DopplerCool
from LAX_exp.system.subsequences.sideband_cool_pulsed import SidebandCoolPulsed
from LAX_exp.system.subsequences.sideband_cool_continuous import SidebandCoolContinuous
from LAX_exp.system.subsequences.sideband_cool_continuous_RAM import SidebandCoolContinuousRAM
__all__.extend(['DopplerCool', 'SidebandCoolPulsed', 'SidebandCoolContinuous', 'SidebandCoolContinuousRAM'])

# state preparation
from LAX_exp.system.subsequences.spin_polarization import SpinPolarization
from LAX_exp.system.subsequences.initialize_qubit import InitializeQubit
__all__.extend(['SpinPolarization', 'InitializeQubit'])

# motional state manipulation
from LAX_exp.system.subsequences.squeeze import Squeeze
from LAX_exp.system.subsequences.squeeze_configurable import SqueezeConfigurable
from LAX_exp.system.subsequences.displace import Displace
from LAX_exp.system.subsequences.tickle_dds import TickleDDS
from LAX_exp.system.subsequences.tickle_fast_dds import TickleFastDDS
from LAX_exp.system.subsequences.tickle_fast_phaser import TickleFastPhaser
from LAX_exp.system.subsequences.parametric_excite import ParametricExcite
from LAX_exp.system.subsequences.QVSA_pulse import QVSAPulse
__all__.extend(['Squeeze', 'SqueezeConfigurable',
                'Displace', 'TickleDDS', 'TickleFastDDS', 'TickleFastPhaser',
                'ParametricExcite', 'QVSAPulse'])

# spin state manipulation
from LAX_exp.system.subsequences.rabi_flop import RabiFlop
from LAX_exp.system.subsequences.ramsey import Ramsey
from LAX_exp.system.subsequences.qubit_pulseshape import QubitPulseShape
from LAX_exp.system.subsequences.qubit_pulseshape_test import QubitPulseShapeTest
from LAX_exp.system.subsequences.qubit_RAP import QubitRAP
__all__.extend(['RabiFlop', 'Ramsey', 'QubitPulseShape', 'QubitRAP', 'QubitPulseShapeTest'])

# readout
from LAX_exp.system.subsequences.readout import Readout
from LAX_exp.system.subsequences.readout_adaptive import ReadoutAdaptive
from LAX_exp.system.subsequences.sideband_readout import SidebandReadout
from LAX_exp.system.subsequences.rabiflop_readout import RabiflopReadout
from LAX_exp.system.subsequences.absorption_probe import AbsorptionProbe
from LAX_exp.system.subsequences.absorption_probe2 import AbsorptionProbe2
from LAX_exp.system.subsequences.doppler_recooling import DopplerRecooling
from LAX_exp.system.subsequences.fock_state_generator import FockStateGenerator
__all__.extend(['Readout', 'ReadoutAdaptive', 'SidebandReadout', 'RabiflopReadout',
                'AbsorptionProbe', 'AbsorptionProbe2', 'DopplerRecooling',
                'FockStateGenerator'])
