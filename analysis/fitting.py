"""
LAX.analysis.fitting

Contains modules used for fitting datasets.
"""

__all__ = ['fitRabiFlopping', 'fitSinc']


# necessary imports
import numpy as np

from scipy.optimize import curve_fit
# todo: separate damped exponential fit from rabi flopping fit


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
    # use damped harmonic oscillator for simplicity
    def fit_func(x, a, b, c):
        """
        todo: document arguments
        """
        return 0.5 * (1. - a * np.exp(-b * x) * np.cos(c * x))

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # todo: these starting params assume data starts at t=0; make agnostic
    # get position of first peak and use as pi/2 time
    max_ind0 =          np.argmax(data_y)
    t0, a0 =            data[max_ind0]
    # get position of next peak and use to find decay constant
    points_margin =     int(np.round(max_ind0 * 1.5))
    max_ind1 =          np.argmax(data_y[points_margin])
    t1, a1 =            data[points_margin + max_ind1]

    # create array of initial guess parameters
    a_guess =           np.abs(1. - 2. * a0)
    b_guess =           0.5 * np.log(a1 / a0) / (t0 - t1)
    c_guess =           np.pi / t0
    param_guess =       np.array([a_guess, b_guess, c_guess])

    # fit and convert covariance matrix to error (1 stdev)
    param_fit, param_cov =  curve_fit(fit_func, data_x, data_y, param_guess)
    param_err =             np.sqrt(np.diag(param_cov))

    return param_fit, param_err


def fitSinc(data, time_fit_s):
    """
    Fit sinc profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # sinc profile fitting function
    def fit_func(x, a, b, c):
        """
        todo: document arguments
        """
        # return (np.power(
        #     (a * np.pi * time_fit_s) * np.sinc(() * np.sqrt(np.power)),
        #     2.) + d)
        return ((a**2. / (a**2. + (x - b)**2.)) * np.sin((np.pi * time_fit_s) * np.sqrt(a**2. + (x - b)**2.))**2. + c)

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # get position of max as line center
    b0, a0 =            data[np.argmax(data_y)]
    # get offset as average of (median, mean) of data
    c0 =                np.mean([np.median(data_y), np.mean(data_y)])
    # create array of initial guess parameters
    a_guess =           (2. / (np.pi * time_fit_s)) * np.arcsin(np.sqrt(a0 - c0))
    param_guess =       np.array([a_guess, b0, c0])

    # fit and convert covariance matrix to error (1 stdev)
    param_fit, param_cov =  curve_fit(fit_func, data_x, data_y, param_guess)
    param_err =             np.sqrt(np.diag(param_cov))

    return param_fit, param_err


def fitSpectralLine(data):
    """
    Fit spectral line.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass
