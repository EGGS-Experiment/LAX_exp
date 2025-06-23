"""
LAX.extensions.utilities

Contains useful functions for experiments.
"""
__all__ = []


from numpy import int64, array
from numpy.random import shuffle

from collections.abc import Iterable
from itertools import product, zip_longest
from scipy.interpolate import Akima1DInterpolator


'''
Dataset interpolation/calibration
'''
__all__.extend(['interpolate_dataset'])

def interpolate_dataset(target_dataset, calibration_dataset):
    """
    todo: document
    """
    # create an interpolation curve to map x to y
    calib_curve = Akima1DInterpolator(calibration_dataset[:, 0], calibration_dataset[:, 1])

    # return the calibrated dataset
    return array(calib_curve(target_dataset))


'''
Experiment configuration tools
'''
__all__.extend(['flatten_maybe_tuple', 'create_experiment_config', 'riffle'])

def flatten_maybe_tuple(xs):
    """
    Use generator to flatten list with SOME tuples.
    """
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten_maybe_tuple(x)
        else:
            yield x

def create_experiment_config(*args, config_type=int64, shuffle_config=True):
    """
    Creates a 2D array of experiment configuration values.
    Each row represents a different experiment configuration to run.
    :param args: parameter lists to merge into an experiment config.
    :param dtype_: data type of
    :param shuffle_config: whether to shuffle the configuration list.
    :return: todo document
    """
    # create an array of values for the experiment to sweep
    exp_config = array([
        list(flatten_maybe_tuple(vals))
        for vals in product(*args)
    ], dtype=config_type)

    if shuffle_config is True: shuffle(exp_config)
    return exp_config

def riffle(*args):
    """
    Riffles an arbitrary number of arbitrary-length iterables sequentially.
    Iterables can each be of different length.
    :param args: iterables to riffle.
    :return: a tuple of the riffled iterables.
    """
    return tuple(
        y
        for x in zip_longest(*args, fillvalue=None)
        for y in x if y is not None
    )

