"""
LAX_exp.system.objects.PulseShaper

Amplitude pulse shaping for waveforms.
For more information, see here: https://en.wikipedia.org/wiki/Window_function.

Note: all of these functions only calculate the rising edge of the given window.
"""
__all__ = ['available_pulse_shapes']

import math
import numpy as np


'''
Square-ish windows
'''
def square(time_arr, time_rolloff, kwargs={}):
    """
    Square pulse.
    For customizability.
    """
    return np.ones(len(time_arr))


'''
Cosine-sum windows
'''
def sine_squared(time_arr, time_rolloff, kwargs={}):
    """
    Sine-squared pulse shape (rising edge only).
    Also known as the Hann or raised-cosine window.
    """
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate sine squared window
    return np.sin(x_vals_readjusted) ** 2.

def blackman(time_arr, time_rolloff, kwargs={}):
    """
    Approximate blackman pulse shape (rising edge only).
    """
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate window
    return (
            0.42 -
            0.5 * np.cos(2. * x_vals_readjusted) +
            0.08 * np.cos(4. * x_vals_readjusted)
    )

def nuttall(time_arr, time_rolloff, kwargs={}):
    """
    Nuttall window (continuous first derivative) (rising edge only).
    """
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate window
    return (
            0.355768 +
            -0.487396 * np.cos(2. * x_vals_readjusted) +
            0.144232 * np.cos(4. * x_vals_readjusted) +
            -0.012604 * np.cos(6. * x_vals_readjusted)
    )

def generalized_cosine(time_arr, time_rolloff, kwargs={}):
    """
    Generalized cosine window (rising edge only).
    kwargs:
        coefficients: a list of coefficients for each term. Should be nonnegative.
    """
    # retrieve coefficients - default is [a0 = 0]
    coeffs = kwargs.get("coefficients", [0])

    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr

    # create window based off coefficients
    # todo: maybe normalize window vals? idk
    return np.sum([
        ((-1)**i) * coeffs[i] * np.cos(2 * i * x_vals_readjusted)
        for i in range(len(coeffs))
    ], axis=0)
    # print(np.shape(np.sum([
    #     ((-1) ** i) * coeffs[i] * np.cos(i * x_vals_readjusted)
    #     for i in range(len(coeffs))
    # ], axis=0)))
    # print(np.sum([
    #     ((-1)**i) * coeffs[i]  * np.cos(i * x_vals_readjusted)
    #     for i in range(len(coeffs))
    # ], axis=1))


'''
Gaussian-type windows
'''
def gaussian(time_arr, time_rolloff, kwargs={}):
    """
    Gaussian pulse shape (rising edge only).
    """
    # get stdev & ensure valid
    std = kwargs.get("std", 0.4)
    if not (0. < std <= 0.5):
        raise ValueError("Error in Gaussian pulse shaping."
                         "stdev must be in (0., 0.5].")
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = 1. / (std * time_rolloff)
    x_vals_readjusted = scale_factor_x * (time_arr - time_rolloff)
    # calculate gaussian window
    return np.exp(-0.5 * (x_vals_readjusted) ** 2.)

def error_function(time_arr, time_rolloff, kwargs={}):
    """
    Error-function pulse shape (rising edge only).
    """
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (2. * np.pi) / time_rolloff
    x_vals_readjusted = (scale_factor_x * time_arr) - np.pi
    # calculate erf window
    return np.array([(1. + math.erf(x_val)) / 2. for x_val in x_vals_readjusted])


'''
Sinc-type windows
'''
def flat_top(time_arr, time_rolloff, kwargs={}):
    """
    Flat top window (rising edge only).
    # todo: actually implement
    """
    # rescale x-axis to the rolloff-time (rising edge only)
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate window
    return (
            0.42 -
            0.5 * np.cos(2. * x_vals_readjusted) +
            0.08 * np.cos(4. * x_vals_readjusted)
    )


available_pulse_shapes = {
    "square": square,

    "sine_squared": sine_squared,
    "blackman": blackman,
    "nuttall": nuttall,
    "generalized_cosine": generalized_cosine,

    "gaussian": gaussian,
    "error_function": error_function,

    "flat_top": flat_top
}


if __name__ == "__main__":
    """
    Test pulse shaping
    # todo: use list of available_pulse_shapes instead of generating our own shape_list lmao
    """
    import matplotlib.pyplot as plt

    # prepare values
    _time_rollon_arb = 100
    _x_vals = np.linspace(0., 200., 100)
    shape_list = [
        # {
        #     "time_rollon_arb": _time_rollon_arb,
        #     "x_vals": _x_vals,
        #     "pulse_shape": 'square'
        # },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "pulse_shape": 'sine_squared'
        },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "pulse_shape": 'blackman'
        },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "pulse_shape": 'nuttall'
        },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "kwargs": {
                "coefficients": [0.42, 0.5, 0.08] # approximate blackman coefficients
            },
            "pulse_shape": 'generalized_cosine'
        },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "pulse_shape": 'error_function'
        },
        {
            "time_rollon_arb": _time_rollon_arb,
            "x_vals": _x_vals,
            "pulse_shape": 'gaussian'
        }
        # {
        #     "time_rollon_arb": _time_rollon_arb,
        #     "x_vals": _x_vals,
        #     "pulse_shape": 'flat_top'
        # }
    ]

    # calculate pulse shapes
    for ps_dict in shape_list:
        ps_dict.update({
            "y_vals": available_pulse_shapes[ps_dict["pulse_shape"]](
                ps_dict["x_vals"],
                ps_dict["time_rollon_arb"],
                ps_dict.get("kwargs", {})
            )
        })

    # plot! - time series
    for ps_dict in shape_list:
        plt.plot(
            ps_dict["x_vals"], ps_dict["y_vals"], '-o', ms=3,
            label="{:s} @ {:4g}".format(ps_dict["pulse_shape"], ps_dict["time_rollon_arb"])
        )
    plt.title("Pulse Shape Testing\nTime Series")
    plt.ylim(0., 1)
    plt.xlabel("Time (arb)")
    plt.ylabel("Amplitude (fractional)")
    plt.legend()
    plt.grid(visible=True)
    plt.show()

    # todo: make it spectral leakage instead of just FFT, since FFT of just the window (no cosine is kinda useless)
    # process FFT
    for ps_dict in shape_list:
        # prepare x-axis for FFT
        x_vals = ps_dict["x_vals"]
        freq_resolution = 1. / (x_vals[-1] - x_vals[0])
        max_freq = (len(x_vals) - 1) * freq_resolution

        # calculate FFT and update dict
        print("orig: {}\tnew: {}".format(len(x_vals), len(np.arange(0., 0.5 * max_freq, freq_resolution))))
        ps_dict.update(
            {
                "fft_x_vals": np.arange(0., 0.5 * max_freq, freq_resolution),
                "fft_y_vals": np.abs(np.fft.fft(ps_dict["y_vals"], norm="forward")[:len(x_vals) // 2]) ** 2.
            }
        )

    # plot! - FFT
    for ps_dict in shape_list:
        plt.plot(
            ps_dict["fft_x_vals"], ps_dict["fft_y_vals"], '-o', ms=3,
            label="{:s} @ {:4g}".format(ps_dict["pulse_shape"], ps_dict["time_rollon_arb"])
        )
    plt.title("Pulse Shape Testing\nFFT")
    plt.xlabel("Frequency (bins)")
    # plt.xscale("log")
    plt.ylabel("Power (dB)")
    plt.yscale("log")
    plt.legend()
    plt.grid(visible=True)
    plt.show()

