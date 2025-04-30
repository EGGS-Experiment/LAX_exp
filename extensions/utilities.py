"""
LAX.extensions.utilities

Contains useful functions for experiments.
# todo: create config creator function
# todo: create interpolator for calibration datasets
# todo: create period rounding multiple function
"""

# __all__ = ["create_experiment_config"]
__all__ = []

# necessary imports
import numpy as np
from itertools import product
from collections.abc import Iterable


# # general
# seconds_to_mu =             lambda seconds:         int64(seconds * 1.e9)
# us_to_mu =                  lambda seconds:         int64(seconds * 1.e3)
# __all__.extend(['seconds_to_mu', 'us_to_mu'])
#
# # Urukul conversions
# hz_to_ftw =                 lambda hz:              int32(round(hz / 1.e9 * 0xFFFFFFFF))
# mhz_to_ftw =                lambda mhz:             int32(round(mhz / 1.e3 * 0xFFFFFFFF))
# pct_to_asf =                lambda pct:             int32(round(pct / 100. * 0x3FFF))
# att_to_mu =                 lambda att:             int32(0xFF) - int32(round(att * 8.))
# __all__.extend(['hz_to_ftw', 'mhz_to_ftw', 'pct_to_asf', 'att_to_mu'])

# def create_experiment_config(*args, shuffle=True):
#     """
#     todo: document
#     Arguments:
#         *args: parameter lists to merge into an experiment config.
#         shuffle: whether to shuffle the configuration list.
#     """
#     def flatten(xs):
#         for x in xs:
#             if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
#                 yield from flatten(x)
#             else:
#                 yield x
#
#     # create an array of values for the experiment to sweep
#     exp_config = np.array([
#         list(flatten(vals))
#         for vals in product(
#             freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
#             vals_pulse4_mu_pow_list
#         )
#     ], dtype=np.int64)
#
#     if shuffle is True:
#         np.random.shuffle(exp_config)
#
#     return exp_config

