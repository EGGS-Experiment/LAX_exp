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
from LAX_exp.system.devices.beam_tickle import BeamTickle
__all__.extend([
    'Beam397Probe', 'Beam397Pump', 'Beam866',
    'Beam854', 'Beam729',
    'BeamTickle'
])

# PMT
from LAX_exp.system.devices.pmt import PMTCounter
__all__.extend(['PMTCounter'])

# TTL triggers
from LAX_exp.system.devices.trigger_linetrigger import Linetrigger
from LAX_exp.system.devices.trigger_rf_modulation import RFSync
__all__.extend(['Linetrigger', 'RFSync'])
