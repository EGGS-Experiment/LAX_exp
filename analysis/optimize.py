"""
LAX.analysis.optimize

Contains modules used for optimization.
"""

__all__ = ['complexLinearFitMinimize', 'complexFitMinimize', 'complexParametricFitMinimize']


# necessary imports
import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares, lsq_linear

# todo: write an optimizer/gradient descent module


'''
Linear Optimization
'''
def complexLinearFitMinimize(dataset):
    """
    Extract the optimal voltage to minimize complex displacement
    from the RF origin.

    Arguments:
        dataset (list(list(float, complex)): the dataset comprised of the real
            independent variable, and the complex dependent variable.

    Returns:
        (float) : the value of the independent variable that minimizes the complex amplitude.
    """
    # create y-vector
    vectorY = dataset[:, 1]
    # create design matrix
    matrixA = np.array([np.ones(len(vectorY)), dataset[:, 0]]).transpose()

    # do a complex LINEAR least squares fit
    res = lsq_linear(matrixA, vectorY)
    b_fit_re, b_fit_im = (res.x[0].real, res.x[0].imag)
    m_fit_re, m_fit_im = (res.x[1].real, res.x[1].imag)
    # print('\t\t\tb_param: {:.3f} + i * {:.3f}'.format(b_fit_re, b_fit_im))
    # print('\t\t\tm_param: {:.3f} + i * {:.3f}\n'.format(m_fit_re, m_fit_im))

    # todo: get error (result.fun; vector of residuals at the soln)
    # extract optimal voltage to minimize displacement
    voltage_optimal = - (b_fit_re * m_fit_re + b_fit_im * m_fit_im) / (m_fit_re**2. + m_fit_im**2.)
    return voltage_optimal


def complexFitMinimize(dataset):
    """
    Extract the optimal voltage to minimize complex displacement
    from the RF origin.

    Arguments:
        dataset (list(list(float, complex)): the dataset comprised of the real
            independent variable, and the complex dependent variable.

    Returns:
        (float) : the value of the independent variable that minimizes the complex amplitude.
    """
    # create norm function for least squares optimization
    def func_norm(b_params, x, y):
        return b_params[0] + b_params[1] * x - y

    # create wrapper function for least_squares to handle complex numbers
    def func_wrap(b_params, x, y):
        # convert tuples to complex numbers
        _b_params = ((b_params[0] + 1.j*b_params[1]),
                     (b_params[2] + 1.j*b_params[3]))
        _y =        y[0] + 1.j*y[1]

        # calculate norm
        res_norm = func_norm(_b_params, x[3], _y)
        return np.array([res_norm.real, res_norm.imag])

    # guess slope as (y_max - y_min) / (x_max - x_min)
    m_guess = (dataset[0, 1] - dataset[-1, 1]) / (dataset[0, 0] - dataset[-1, 0])
    # guess x-intercept as mean of y - mx
    b_guess = np.median(dataset[:, 1] - m_guess * dataset[:, 0])
    params_guess = (b_guess.real, b_guess.imag, m_guess.real, m_guess.imag)

    # debug printouts
    print('\t\tguess b_param: {:.3f} + i * {:.3f}'.format(b_guess.real, b_guess.imag))
    print('\t\tguess m_param: {:.3f} + i * {:.3f}'.format(m_guess.real, m_guess.imag))

    # do a complex least squares fit
    res = least_squares(func_wrap, params_guess, args=(dataset[:, 0], dataset[:, 1]), jac='3-point', method='lm')
    b_fit_re, b_fit_im, m_fit_re, m_fit_im = res.x
    print('\t\t\tb_param: {:.3f} + i * {:.3f}'.format(b_fit_re, b_fit_im))
    print('\t\t\tm_param: {:.3f} + i * {:.3f}\n'.format(m_fit_re, m_fit_im))

    # todo: get error

    # extract optimal voltage to minimize displacement
    voltage_optimal = - (b_fit_re * m_fit_re + b_fit_im * m_fit_im) / (m_fit_re**2. + m_fit_im**2.)
    return voltage_optimal


def complexParametricFitMinimize(dataset):
    """
    Extract the optimal voltage to minimize complex displacement
    from the RF origin.
    # todo: note this is deprecated

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
    fit_complex = linregress(dataset_y[:, 0], dataset_y[:, 1])
    m_c, b_c = fit_complex.slope, fit_complex.intercept

    # get the projection of the DV onto the fitted line
    vec_complex = np.array([1, m_c]) / np.linalg.norm(np.array([1, m_c]))
    dataset_proj = np.dot(dataset_y, vec_complex)

    # fit the x dataset to the parameterized line
    fit_parameterized = linregress(dataset_x, dataset_proj)
    m_p, b_p = fit_parameterized.slope, fit_parameterized.intercept

    # extract x_min
    x_min = - (m_c * b_c) / (1 + np.power(m_c, 2))

    # convert x_min to V_min
    V_min = (x_min - b_p) / m_p
    return V_min
