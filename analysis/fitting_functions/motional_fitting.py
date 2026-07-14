from numpy import interp, linspace, exp

from scipy.signal import find_peaks
from scipy.special import laguerre
from scipy.optimize import curve_fit

def phonon_conversion(alpha):
    return 1 - exp(-alpha**2)*(laguerre(0)(alpha**2))**2

def nbar(alpha):
    return alpha**2

def fit_rap(state_vals_ave):

    alphas = linspace(0, 5,100)
    nbars = [nbar(alpha) for alpha in alphas]

    phonon_conversions = [phonon_conversion(alpha) for alpha in alphas]
    phonons = interp(state_vals_ave, phonon_conversions, nbars)
    return phonons