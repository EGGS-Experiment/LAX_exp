"""
LAX.extensions.constants

Important numerical motifs for hardware, software, and experiment.
"""

__all__ = ['TIME_PROFILESWITCH_DELAY_MU', 'TIME_PHASEAUTOCLEAR_DELAY_MU', 'TIME_RFSWITCH_DELAY_MU']


# urukul AD9910 timings
# todo: change profileswitch timings to reflect actual delay between end of cfg_write and actual update; should be ~79ns just like for phase autoclear
TIME_PROFILESWITCH_DELAY_MU =   420
TIME_PHASEAUTOCLEAR_DELAY_MU =  78
TIME_RFSWITCH_DELAY_MU =        31

# todo: add ARTIQ phaser timings as well
