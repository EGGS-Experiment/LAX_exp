"""
Contains base device classes for LAX.
"""

__all__ = []

# beams
from .beam_urukul import Beam_Urukul
__all__.extend(["Beam_Urukul"])

# TTLs
from .trigger_ttl import Trigger_TTL
__all__.extend(["Trigger_TTL"])
