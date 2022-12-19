"""
LAX.extensions.conversions

Contains useful conversion functions.
"""

__all__ = []


# necessary imports
from numpy import int32, int64

# general
seconds_to_mu =             lambda seconds:         int64(seconds * 1e9)
us_to_mu =                  lambda seconds:         int64(seconds * 1e3)
__all__.extend(['seconds_to_mu', 'us_to_mu'])

# DDS parameters
mhz_to_ftw =                lambda mhz:             int32(round(mhz / 1e3 * 0xFFFFFFFF))
pct_to_asf =                lambda pct:             int32(round(pct / 100 * 0x3FFF))
__all__.extend(['mhz_to_ftw', 'pct_to_asf'])
