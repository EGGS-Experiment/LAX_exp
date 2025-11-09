"""
LAX.extensions.constants

Important numerical motifs for hardware, software, and experiment.
"""

__all__ = []


from numpy import int32, int64


'''
DATA TYPES & FORMATS
'''
TYPE_DMA_HANDLE = (0, int64(0), int32(0), False)
__all__.extend(['TYPE_DMA_HANDLE'])


'''
TIMINGS
'''
# AD9910 timings
TIME_AD9910_PROFILE_SWITCH_DELAY_MU =   63
TIME_AD9910_PHASE_AUTOCLEAR_DELAY_MU =  79
# todo: set_mu total delay
__all__.extend(['TIME_AD9910_PROFILE_SWITCH_DELAY_MU', 'TIME_AD9910_PHASE_AUTOCLEAR_DELAY_MU'])

# Urukul timings
TIME_URUKUL_BUS_WRITE_DELAY_MU =    416
TIME_URUKUL_RFSWITCH_DELAY_MU =     140
# todo: att write delay
__all__.extend(['TIME_URUKUL_BUS_WRITE_DELAY_MU', 'TIME_URUKUL_RFSWITCH_DELAY_MU'])

# Phaser timings
TIME_PHASER_OSCILLATOR_PIPELINE_LATENCY_MU =    1953
__all__.extend(['TIME_PHASER_OSCILLATOR_PIPELINE_LATENCY_MU'])

# external device timings
TIME_ZASWA2_SWITCH_DELAY_MU =   31
__all__.extend(['TIME_ZASWA2_SWITCH_DELAY_MU'])

