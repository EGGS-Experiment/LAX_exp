"""
LAX.analysis.fitting

Contains modules used for fitting datasets.
"""
__all__ = [
    'fitDampedOscillator', 'fitDampedDrivenOscillatorAmplitude', 'fitDampedDrivenOscillatorPhase',
    'fitSinc', 'fitSincGeneric', 'fitGaussian', 'fitLorentzian', 'fitVoigt',
    'fitLine', 'fitLineLinear'
]

import numpy as np
from scipy.optimize import curve_fit, least_squares, lsq_linear


'''
Fitting: Oscillators
'''
class fitDampedOscillator:
    """
    Fit exponentially damped harmonic oscillator to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # use damped harmonic oscillator for simplicity
    def fit_func(self, x, a, b, c):
        """
        todo: document arguments
        """
        return 0.5 * (1. - a * np.exp(-b * x) * np.cos(c * x))

    def fit(self, data):
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
        if points_margin < len(data):
            max_ind1 = np.argmax(data_y[points_margin])
            t1, a1 = data[points_margin + max_ind1]
        else:
            t1 = t0 * 3
            a1 = a0

        # create array of initial guess parameters
        a_guess = np.abs(1. - 2. * a0)
        b_guess = 0.5 * np.log(a1 / a0) / (t0 - t1)
        c_guess = np.pi / t0
        param_guess = np.array([a_guess, b_guess, c_guess])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))
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


class fitDampedDrivenOscillatorAmplitude:
    """
    Fit amplitude response of the damped driven harmonic oscillator.

    Note: amplitude response lineshape of the damped driven harmonic oscillator
    is also called a Breit-Wigner profile.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # Breit-Wigner profile
    def fit_func(self, x, a, b, c):
        """
        todo: document arguments
        """
        return a / ((x ** 2. - b ** 2.) ** 2. + (x * c) ** 2.) ** 0.5

    def fit(self, data):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # extract starting parameter guesses
        # guess linewidth as 2 kHz
        c_guess = 0.002
        # get position and value of peak
        b0, a0 = data[np.argmax(data_y)]
        # since peak doesn't occur at linecenter for breit-wigner,
        # convert peak values to actual guess values
        a_guess = 0.5 * a0 * c_guess
        b_guess = (2 * b0 ** 2. - c_guess) ** 0.5
        # create array of initial guess parameters
        param_guess = np.array([a_guess, b_guess, c_guess])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))
        return param_fit, param_err


class fitDampedDrivenOscillatorPhase:
    """
    Fit phase response of the damped driven harmonic oscillator.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # regular old lorentzian
    def fit_func(self, x, a, b, mu):
        """
        todo: document arguments
        """
        return a / ((x - mu) ** 2. + (b) ** 2.)

    def fit(self, data):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # extract starting parameter guesses
        # get position of max as linecenter
        mu0, a0 = data[np.argmax(data_y)]
        # get b0 by numerically guessing FWHM
        b0 = np.abs(data_x[np.argmin(np.abs(data_y - 0.5 * a0))] - mu0)
        # create array of initial guess parameters
        param_guess = np.array([a0, b0, mu0])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))

        # convert b_fit to \Gamma / 2 (i.e. linewidth)
        param_fit[1] *= 2.
        param_err[1] *= 2.
        return param_fit, param_err


'''
Fitting: Spectroscopy Profiles
'''


class fitSinc:
    """
    Fit Sinc profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    def __init__(self):
        self.time_fit_s = None

    # use sinc profile with lorentzian amplitude per optical bloch equations
    def fit_func(self, x, a, b, c):
        """
        todo: document arguments
        """
        return ((a ** 2. / (a ** 2. + (x - b) ** 2.)) * np.sin(
            (np.pi * self.time_fit_s) * (a ** 2. + (x - b) ** 2.) ** 0.5) ** 2. + c)

    def fit(self, data, time_fit_s):

        # assign time_s
        self.time_fit_s = time_fit_s

        # separate data into x and y
        try:
            data = np.array(data)
            data_x, data_y = data.transpose()
        except:
            data = np.array(data)
            data_x = data[0, :]
            data_y = data[1, :]

        # extract starting parameter guesses
        b0 = np.mean(data_x[np.argwhere(data_y == np.max(data_y))])
        # get offset as average of (median, min) of data
        c0 = 0.5 * (np.median(data_y) + np.min(data_y))
        # guess amplitude using max y-value with offset subtracted
        a0 = np.max(data_y) - c0
        a0 = ((2. / (np.pi * time_fit_s)) * np.arcsin(np.sqrt(a0)))

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, p0=[a0, b0, c0])
        param_err = np.sqrt(np.diag(param_cov))
        return param_fit, param_err


class fitGaussian:
    """
    Fit gaussian profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # regular old gaussian
    def fit_func(self, x, a, b, mu):
        """
        todo: document arguments
        """
        return a * np.exp(-b * (x - mu) ** 2.)

    def fit(self, data):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # extract starting parameter guesses
        # get position of max as linecenter
        mu0, a0 = data[np.argmax(data_y)]
        # get b0 by numerically guessing FWHM
        gamma0 = np.abs(data_x[np.argmin(np.abs(data_y - 0.5 * a0))] - mu0)
        b0 = np.log(2) / gamma0 ** 2.
        # create array of initial guess parameters
        param_guess = np.array([a0, b0, mu0])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))
        return param_fit, param_err


class fitLorentzian:
    """
    Fit Lorentzian profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # regular old lorentzian
    def fit_func(self, x, a, b, mu):
        """
        todo: document arguments
        """
        return a / ((x - mu) ** 2. + (b) ** 2.)

    def fit(self, data):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # extract starting parameter guesses
        # get position of max as linecenter
        mu0, a0 = data[np.argmax(data_y)]
        # get b0 by numerically guessing FWHM
        b0 = np.abs(data_x[np.argmin(np.abs(data_y - 0.5 * a0))] - mu0)
        # create array of initial guess parameters
        param_guess = np.array([a0, b0, mu0])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))

        # convert b_fit to \Gamma / 2 (i.e. linewidth)
        param_fit[1] *= 2.
        param_err[1] *= 2.
        return param_fit, param_err


class fitVoigt:
    """
    Fit Voigt profile to data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # voigt profile: convolution of lorentzian (e.g. spectrum) with gaussian (e.g. temperature)
    def fit_func(self, x, a, b, mu):
        """
        todo: document arguments
        """
        return a * np.exp(-b * (x - mu) ** 2.)

    def fit(self, data):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # extract starting parameter guesses
        # get position of max as line center
        mu0, a0 = data[np.argmax(data_y)]
        # get b0 by numerically guessing FWHM
        gamma0 = np.abs(data_x[np.argmin(np.abs(data_y - 0.5 * a0))] - mu0)
        b0 = np.log(2) / gamma0 ** 2.
        # create array of initial guess parameters
        param_guess = np.array([a0, b0, mu0])

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, data_x, data_y, param_guess)
        param_err = np.sqrt(np.diag(param_cov))
        return param_fit, param_err


'''
Fitting: Other
'''


class fitLine:
    """
    Fit linear trend to data using nonlinear least-squares.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    # create norm function for least squares optimization
    def func_norm(self, b_params, x, y):
        """
        todo: document arguments
        """
        # todo: shouldn't this be different? why is it x[1]?
        # todo: welp nope; yes it should be x[1]; least_squares passes in an arr matching b_params dimension
        return b_params[0] + (b_params[1] * x[1]) - y

    def fit(self, data, bounds=(-np.inf, np.inf)):
        # guess starting parameters for the line fit
        # guess slope as (y_max - y_min) / (x_max - x_min)
        m_guess = (data[0, 1] - data[-1, 1]) / (data[0, 0] - data[-1, 0])
        # guess x-intercept as median of y - mx
        b_guess = np.median(data[:, 1] - m_guess * data[:, 0])

        # do a linear least squares fit to the data
        res = least_squares(self.func_norm, [b_guess, m_guess], args=(data[:, 0], data[:, 1]), bounds=bounds)
        res_intercept, res_slope = res.x
        # todo: make it follow param_fit, param_err return style
        return res_intercept, res_slope


class fitLineLinear:
    """
    Fit linear trend to data using linear least-squares.

    Arguments:
        ***todo

    Returns:
        ***todo
    """

    def fit(self, data, bounds=(-np.inf, np.inf)):
        # separate data into x and y
        data = np.array(data)
        data_x, data_y = data.transpose()

        # create design matrix
        matrixA = np.array([np.ones(len(data_y)), data_x]).transpose()

        # do a linear least squares fit
        res = lsq_linear(matrixA, data_y)
        param_fit = res.x
        return param_fit


class fitSincGeneric:
    """
    Fitting function for a generic sinc function

    Args:
        x: x values
        y: y values

    Returns:
        param_fit: optimized parameters from curve_fitting
        param_err: estimated covariance from curve_fitting
        ydata: generated curve over given domain
    """

    def fit_func(self, x, a, b, c, d):
        """
        Fitting function for a generic sinc function

        Args:
            a: amplitude
            b: linecenter
            c: linewidth
            d: amplitude offset

        Returns:
            sinc function
        """
        return a * np.sinc(c * (x - b)) ** 2. + d

    def fit(self, x: np.array, y: np.array, amplitude_factor=1):
        ## extract starting parameter guesses
        # get indices of max y-values
        indices_max_y = np.argwhere(y == np.max(y))
        # get offset as average of (median, min) of data
        b0 = np.mean(x[indices_max_y])
        # get offset as min of data
        d0 = np.min(y)
        # guess amplitude using max y-value with offset subtracted
        a0 = amplitude_factor * (np.max(y) - d0)
        # numerically guess FWHM for a sinc function
        c0 = 0.69578 / np.abs(x[np.argmin(np.abs(y - 0.5 * np.max(y)))] - b0)

        # fit and convert covariance matrix to error (1 stdev)
        param_fit, param_cov = curve_fit(self.fit_func, x, y, p0=[a0, b0, c0, d0])
        xdata = np.linspace(np.min(x), np.max(x), int(1e6))
        ydata = self.fit_func(xdata, *param_fit)
        param_err = np.sqrt(np.diag(param_cov))
        return param_fit, param_err, ydata
