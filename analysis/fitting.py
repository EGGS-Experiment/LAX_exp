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
        return 0.5 * (1. - a * np.exp(-b * x) * np.cos(c * x))

    # separate data into x and y
    data = np.array(data)
    data_x, data_y = data.transpose()

    # extract starting parameter guesses
    # todo: these starting params assume data starts at t=0; make agnostic
    # get position of first peak and use as pi/2 time
    max_ind0 = np.argmax(data_y)
    t0, a0 = data[max_ind0]
    # get position of next peak and use to find decay constant
    points_margin = int(np.round(max_ind0 * 1.5))
    max_ind1 = np.argmax(data_y[points_margin])
    t1, a1 = data[points_margin + max_ind1]

    # create array of initial guess parameters
    a_guess = np.abs(1. - 2. * a0)
    b_guess = 0.5 * np.log(a1 / a0) / (t0 - t1)
    c_guess = np.pi / t0
    param_guess = np.array([a_guess, b_guess, c_guess])
    # print('\tguess: {}'.format(param_guess))

    # fit!
    param_fit, param_cov = curve_fit(fit_func, data_x, data_y, param_guess)
    # convert coveriance matrix to error (1 stdev)
    param_err = np.sqrt(np.diag(param_cov))
    # print('\tactual: {}'.format(param_fit))

    return param_fit, param_err


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
