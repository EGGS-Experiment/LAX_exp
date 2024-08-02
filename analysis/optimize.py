"""
LAX.analysis.optimize

Contains modules used for optimization.
"""

__all__ = ['complexLinearFitMinimize']


# necessary imports
import numpy as np
from scipy.optimize import lsq_linear

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
        tuple(float, float): the value and error of the independent variable that minimizes the complex amplitude.
    """
    # create y-vector
    vectorY = dataset[:, 1]
    # create design matrix
    matrixA = np.array([np.ones(len(vectorY)), dataset[:, 0]]).transpose()

    # do a complex LINEAR least squares fit
    res = lsq_linear(matrixA, vectorY)
    b_fit_re, b_fit_im = (res.x[0].real, res.x[0].imag)
    m_fit_re, m_fit_im = (res.x[1].real, res.x[1].imag)

    # extract optimal voltage to minimize displacement
    voltage_optimal =   - (b_fit_re * m_fit_re + b_fit_im * m_fit_im) / (m_fit_re**2. + m_fit_im**2.)
    voltage_err =       np.sum(np.abs(res.fun)) / np.sqrt(len(res.fun))
    # todo: calculate "amplitude" so we can scale in future
    return voltage_optimal, voltage_err

