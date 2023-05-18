"""
LAX.analysis.signals

Contains modules used for signal analysis (e.g. demodulation)
"""

__all__ = []


# necessary imports
import numpy as np
from scipy import stats


'''
Modulation/Demodulation
'''
# @rpc(flags={"async"})
# def demodulateTimestamps(self, freq_hz: TFloat, timestamp_list_s: TArray(TInt32)) -> TList(TFloat):
#     """
#     Demodulate a set of timestamps with respect to some known frequency.
#
#     Arguments:
#         freq_hz     (int)   : the frequency to demodulate at (in Hz)
#     """
#     # remove starting time and digitally demodulate counts
#     counts_mu = self.core.mu_to_seconds(np.array(timestamp_mu_list) - time_start_mu)
#     counts_demod = np.sum(np.exp((2.j * np.pi * freq_mhz * 1e6) * counts_mu)) / self.num_counts
#
#     return counts_demod
