"""
Contains base device classes for LAX.
"""

__all__ = []

# device classes
from .base_device import LAXDevice
__all__.extend(['LAXDevice'])

# tmp remove
from numpy import int32, int64
mhz_to_ftw = lambda mhz: int32(round(mhz / 1e3 * 0xFFFFFFFF))
pct_to_asf = lambda pct: int32(round(pct / 100 * 0x3FFF))
seconds_to_mu = lambda seconds: int64(seconds * 1e9)
__all__.extend(['mhz_to_ftw', 'pct_to_asf'])
