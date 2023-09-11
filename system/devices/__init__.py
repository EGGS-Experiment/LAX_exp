"""
Contains specifically constructed device classes for LAX.
"""

__all__ = []

# beams
from LAX_exp.system.devices.beam_397_probe import Beam397Probe
from LAX_exp.system.devices.beam_397_pump import Beam397Pump
from LAX_exp.system.devices.beam_866 import Beam866
from LAX_exp.system.devices.beam_854 import Beam854
from LAX_exp.system.devices.beam_729 import Beam729
__all__.extend([
    'Beam397Probe', 'Beam397Pump', 'Beam866',
    'Beam854', 'Beam729',
])

# DDS
from LAX_exp.system.devices.dds_modulation import DDSModulation
from LAX_exp.system.devices.dds_dipole import DDSDipole
__all__.extend(['DDSModulation', 'DDSDipole'])

# AWG
from LAX_exp.system.devices.phaser_eggs import PhaserEGGS
__all__.extend(['PhaserEGGS'])

# PMT
from LAX_exp.system.devices.pmt import PMTCounter
__all__.extend(['PMTCounter'])

# TTL triggers
from LAX_exp.system.devices.trigger_line import TriggerLine
from LAX_exp.system.devices.trigger_rf import TriggerRF
__all__.extend(['TriggerLine', 'TriggerRF'])
