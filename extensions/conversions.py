"""
LAX.extensions.conversions

Contains useful conversion functions.
"""

__all__ = []


# necessary imports
from numpy import int32, int64

# general
seconds_to_mu =             lambda seconds:         int64(seconds * 1.e9)
us_to_mu =                  lambda seconds:         int64(seconds * 1.e3)
__all__.extend(['seconds_to_mu', 'us_to_mu'])

# Urukul conversions
hz_to_ftw =                 lambda hz:              int32(round(hz / 1.e9 * 0xFFFFFFFF))
mhz_to_ftw =                lambda mhz:             int32(round(mhz / 1.e3 * 0xFFFFFFFF))
pct_to_asf =                lambda pct:             int32(round(pct / 100. * 0x3FFF))
att_to_mu =                 lambda att:             int32(0xFF) - int32(round(att * 8.))
__all__.extend(['hz_to_ftw', 'mhz_to_ftw', 'pct_to_asf', 'att_to_mu'])

