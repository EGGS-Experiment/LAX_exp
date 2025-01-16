"""
LAX_exp.system.objects.PulseShaper

Amplitude pulse shaping for waveforms.
"""
__all__ = ['available_pulse_shapes']


import math
import numpy as np


def square(time_arr, time_rolloff, kwargs={}):
    """
    Square pulse.
    For customizability.
    """
    return np.ones(len(time_arr))

def sine_squared(time_arr, time_rolloff, kwargs={}):
    """
    Sine-squared pulse shape (rising edge only).
    Also known as the Hann or raised-cosine windows.
    """
    # rescale x-axis to halfway
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate sine squared window
    return np.sin(x_vals_readjusted) ** 2.

def blackman(time_arr, time_rolloff, *wargs):
    """
    Approximate blackman pulse shape (rising edge only).
    """
    # rescale x-axis to halfway
    scale_factor_x = (np.pi / 2.) / time_rolloff
    x_vals_readjusted = scale_factor_x * time_arr
    # calculate window
    return (
            0.42 -
            0.5 * np.cos(2. * x_vals_readjusted) +
            0.08 * np.cos(4. * x_vals_readjusted)
    )

def gaussian(time_arr, time_rolloff, kwargs={}):
    """
    Gaussian pulse shape (rising edge only).
    """
    # get stdev & ensure valid
    std = kwargs.get("std", 0.4)
    if not (0. < std <= 0.5):
        raise ValueError("Error in Gaussian pulse shaping."
                         "stdev must be in (0., 0.5].")
    # rescale x-axis to halfway
    scale_factor_x = 1. / (std * time_rolloff)
    x_vals_readjusted = scale_factor_x * (time_arr - time_rolloff)
    # calculate gaussian window
    return np.exp(-0.5 * (x_vals_readjusted) ** 2.)

def error_function(time_arr, time_rolloff, kwargs={}):
    """
    Error-function pulse shape (rising edge only).
    """
    # rescale x-axis to halfway
    scale_factor_x = (2. * np.pi) / time_rolloff
    x_vals_readjusted = (scale_factor_x * time_arr) - np.pi
    # calculate erf window
    return np.array([(1. + math.erf(x_val)) / 2. for x_val in x_vals_readjusted])


available_pulse_shapes = {
    "square": square,
    "sine_squared": sine_squared,
    "blackman": blackman,
    "gaussian": gaussian,
    "error_function": error_function,
}


if __name__ == "__main__":
    """
    Test pulse shaping
    """
    # prepare values
    shape_list = [
        {
            "time_rollon_arb": 100,
            "x_vals": np.linspace(0., 100., 1000),
            "pulse_shape": 'sine_squared'
        },
        {
            "time_rollon_arb": 100,
            "x_vals": np.linspace(0., 100., 1000),
            "pulse_shape": 'blackman'
        },
        {
            "time_rollon_arb": 100,
            "x_vals": np.linspace(0., 100., 1000),
            "pulse_shape": 'error_function'
        },
        {
            "time_rollon_arb": 100,
            "x_vals": np.linspace(0., 100., 1000),
            "pulse_shape": 'gaussian'
        }
    ]

    # calculate pulse shapes
    for ps_dict in shape_list:
        ps_dict.update({
            "y_vals": available_pulse_shapes[ps_dict["pulse_shape"]](
                ps_dict["x_vals"],
                ps_dict["time_rollon_arb"],
                ps_dict.get("ps_args", {})
            )
        })

    # todo: process FFT

    # plot
    import matplotlib.pyplot as plt
    for ps_dict in shape_list:
        plt.plot(
            ps_dict["x_vals"], ps_dict["y_vals"], '-o', ms=3,
            label="{:s} @ {:4g}".format(ps_dict["pulse_shape"], ps_dict["time_rollon_arb"])
        )
    # tmp remove
    # plt.plot(
    #     shape_list[0]["x_vals"], shape_list[0]["y_vals"]-shape_list[1]["x_vals"], '-o', ms=5,
    #     label="yzde"
    # )
    # tmp remove

    plt.title("Pulse Shape Testing")
    plt.ylim(0., 1)
    # plt.xlim(0., np.inf)
    plt.xlabel("Time (arb)")
    plt.ylabel("Amplitude (fractional)")
    plt.legend()
    plt.grid(visible=True)
    plt.show()

