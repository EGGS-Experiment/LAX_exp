"""
LAX.analysis.signals

Contains modules used for signal analysis (e.g. demodulation)
"""

# __all__ = ["demodulateTimestamps", "complexFitMinimize"]
__all__ = ["complexFitMinimize"]


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

def complexFitMinimize(self, dataset):
    """
    Extract the optimal voltage to minimize complex displacement
    from the RF origin.

    Arguments:
        dataset (list(list(float, complex)): the dataset comprised of the real
            independent variable, and the complex dependent variable.

    Returns:
        (float) : the value of the independent variable that minimizes the complex amplitude.
    """
    # split dataset into IV/DV, and real/imaginary
    dataset_x = dataset[:, 0]
    dataset_y = np.array([np.real(dataset[:, 1]), np.imag(dataset[:, 1])]).transpose()

    # fit the DV in the complex plane and get the unit vector of the line
    fit_complex = stats.linregress(dataset_y[:, 0], dataset_y[:, 1])
    m_c, b_c = fit_complex.slope, fit_complex.intercept

    # get the projection of the DV onto the fitted line
    vec_complex = np.array([1, m_c]) / np.linalg.norm(np.array([1, m_c]))
    dataset_proj = np.dot(dataset_y, vec_complex)

    # fit the x dataset to the parameterized line
    fit_parameterized = stats.linregress(dataset_x, dataset_proj)
    m_p, b_p = fit_parameterized.slope, fit_parameterized.intercept

    # extract x_min
    x_min = - (m_c * b_c) / (1 + np.power(m_c, 2))

    # convert x_min to V_min
    V_min = (x_min - b_p) / m_p

    return V_min
