"""
LAX.analysis.fitting

Contains modules used for fitting datasets.
"""

__all__ = ['fitRabiFlopping']


# necessary imports
import numpy as np

from scipy.optimize import curve_fit


'''
Fitting
'''
def fitRabiFlopping(data):
    """
    Fit Rabi Flopping population data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # fit damped harmonic oscillator for simplicity
    def fit_func(x, a, b, c):
        """
        todo: document arguments
        """
        return a * np.exp(-b * x) * np.power(np.sin(c * x), 2.)

    # separate data into x and y
    data = np.array(data)
    data_x, data_y = data.transpose()

    # extract starting parameter guesses
    # todo: these starting params assume data starts at t=0; make agnostic
    # get position of first peak and use as pi/2 time
    max_ind0 = np.argmax(data_y)
    t0, a0 = data[max_ind0]
    # get position of next peak and use to find decay constant
    points_margin = max_ind0 + np.round(max_ind0 * 2.)
    max_ind1 = np.argmax(data_y[points_margin])
    t1, a1 = data[points_margin + max_ind1]

    # create array of initial guess parameters
    a_guess = a0
    b_guess = np.log(a1 / a0) / (t0 - t1)
    c_guess = np.pi / (2. * t0)
    param_guess = np.array([a_guess, b_guess, c_guess])

    # fit!
    param_fit, param_cov = curve_fit(fit_func, data_x, data_y, param_guess)
    return param_fit


def fitSidebandCooling(data):
    """
    Fit sideband cooling spectrum data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass


def fitSpectralLine(data):
    """
    Fit spectral line.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass
