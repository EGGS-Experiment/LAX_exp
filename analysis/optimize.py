"""
LAX.analysis.optimize

Contains modules used for optimization.
"""

__all__ = ['complexFitMinimize', 'complexParametricFitMinimize']


# necessary imports
import numpy as np
from scipy import stats
from scipy import optimize


'''
Linear Optimization
'''
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
        return b_params[0] + b_params[1] * x[1] - y

    # guess starting b_params
    # todo: is this the best guess we can do? try to improve somehow
    b_guess = [0.1, 0.1]

    # do a complex least squares fit
    res = optimize.least_squares(func_norm, b_guess, args=(dataset[:, 0], dataset[:, 1]))
    res_intercept, res_slope = res.x

    # extract optimal voltage to minimize displacement
    voltage_optimal = - (np.re(res_intercept) * np.re(res_slope) + np.imag(res_intercept) * np.imag(res_slope)) / np.abs(res_slope)**2.

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
