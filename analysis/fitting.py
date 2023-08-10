"""
LAX.analysis.fitting

Contains modules used for fitting datasets.
"""

__all__ = ['fitDampedOscillator', 'fitSinc', 'fitGaussian', 'fitLorentzian', 'fitVoigt', 'fitLine']


# necessary imports
import numpy as np
from scipy.special import wofz
from scipy.optimize import curve_fit, least_squares


'''
Fitting: Oscillators
'''
def fitDampedOscillator(data):
    """
    Fit exponentially damped harmonic oscillator to data.

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


def fitRabiFlopping(data):
    """
    todo: document

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass


'''
Fitting: Spectroscopy Profiles
'''
def fitSinc(data, time_fit_s):
    """
    Fit Sinc profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # use sinc profile with lorentzian amplitude per optical bloch equations
    def fit_func(x, a, b, c):
        """
        todo: document arguments
        """
        return ((a**2. / (a**2. + (x - b)**2.)) * np.sin((np.pi * time_fit_s) * (a**2. + (x - b)**2.)**0.5)**2. + c)

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # get position of max as linecenter
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


def fitGaussian(data):
    """
    Fit gaussian profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # regular old gaussian
    def fit_func(x, a, b, mu):
        """
        todo: document arguments
        """
        return a * np.exp(-b * (x - mu)**2.)

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # get position of max as linecenter
    mu0, a0 =           data[np.argmax(data_y)]
    # get b0 by numerically guessing FWHM
    gamma0 =            data_x[np.argmax(np.abs(data_y - 0.5 * a0))]
    b0 =                np.log(2) / (gamma0 - mu0)**2.
    # create array of initial guess parameters
    param_guess =       np.array([a0, b0, mu0])

    # fit and convert covariance matrix to error (1 stdev)
    param_fit, param_cov =  curve_fit(fit_func, data_x, data_y, param_guess)
    param_err =             np.sqrt(np.diag(param_cov))
    return param_fit, param_err


def fitLorentzian(data):
    """
    Fit Lorentzian profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # regular old lorentzin
    def fit_func(x, a, b, mu):
        """
        todo: document arguments
        """
        return a / ((x - mu)**2. + (b)**2.)

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # get position of max as linecenter
    mu0, a0 =           data[np.argmax(data_y)]
    # get b0 by numerically guessing FWHM
    b0 =                data_x[np.argmax(np.abs(data_y - 0.5 * a0))]
    # create array of initial guess parameters
    param_guess =       np.array([a0, b0, mu0])

    # fit and convert covariance matrix to error (1 stdev)
    param_fit, param_cov =  curve_fit(fit_func, data_x, data_y, param_guess)
    param_err =             np.sqrt(np.diag(param_cov))

    # convert b_fit to \Gamma / 2 (i.e. linewidth)
    param_fit[1] *= 2.
    param_err[1] *= 2.
    return param_fit, param_err


def fitVoigt(data):
    """
    Fit Voigt profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # voigt profile: convolution of lorentzian (e.g. spectrum) with gaussian (e.g. temperature)
    def fit_func(x, a, b, mu):
        """
        todo: document arguments
        """
        return a * np.exp(-b * (x - mu)**2.)

    # separate data into x and y
    data =              np.array(data)
    data_x, data_y =    data.transpose()

    # extract starting parameter guesses
    # get position of max as line center
    mu0, a0 =           data[np.argmax(data_y)]
    # get b0 by numerically guessing FWHM
    gamma0 =            data_x[np.argmax(np.abs(data_y - 0.5 * a0))]
    b0 =                np.log(2) / (gamma0 - mu0)**2.
    # create array of initial guess parameters
    param_guess =       np.array([a0, b0, mu0])

    # fit and convert covariance matrix to error (1 stdev)
    param_fit, param_cov =  curve_fit(fit_func, data_x, data_y, param_guess)
    param_err =             np.sqrt(np.diag(param_cov))
    return param_fit, param_err


'''
Fitting: Other
'''
def fitLine(data):
    """
    Fit linear trend to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # create norm function for least squares optimization
    def func_norm(b_params, x, y):
        """
        todo: document arguments
        """
        return b_params[0] + (b_params[1] * x[1]) - y

    # guess starting parameters for the line fit
    b_guess =           np.min(data[:, 1])
    m_guess =           np.mean(data[:, 1] / data[:, 0])

    # do a linear least squares fit to the data
    res =               least_squares(func_norm, [b_guess, m_guess], args=(data[:, 0], data[:, 1]))
    res_intercept, res_slope = res.x
    return res_intercept, res_slope
