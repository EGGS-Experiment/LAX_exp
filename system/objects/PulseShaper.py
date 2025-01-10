"""
LAX_exp.system.objects.PulseShaper

Amplitude pulse shaping for waveforms.
"""
__all__ = ['available_pulse_shapes']


# necessary imports
import numpy as np
from math import erf

def ps_sine_squared(time_arr, time_rolloff):
    """
    Sine-squared pulse shape (rising edge only).
    """
    # rescale x-axis
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate sine squared window
    return np.sin(x_vals_readjusted) ** 2.

def ps_error_function(time_arr, time_rolloff):
    """
    Error-function pulse shape (rising edge only).
    """
    # rescale x-axis
    scale_factor_x = (2. * np.pi) / time_rolloff
    x_vals_readjusted = (scale_factor_x * time_arr) - np.pi
    # calculate erf window
    return np.array([(1. + erf(x_val)) / 2. for x_val in x_vals_readjusted])

available_pulse_shapes = {
    "sine_squared":     ps_sine_squared,
    "error_function":   ps_error_function
}
