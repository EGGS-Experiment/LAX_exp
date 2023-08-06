"""
LAX.analysis.processing

Contains helpful/commonly used modules for processing datasets.
"""

__all__ = ['findThreshold', 'groupBy']


# necessary imports
import numpy as np

from itertools import groupby

from scipy.stats import iqr
from scipy.signal import find_peaks


'''
Dataset Processing
'''
def findThreshold(counts_arr):
    """
    Get the binary discrimination threshold for a dataset.
    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # calculate histogram bin width using freedman diaconis rule (bin width = 2 * iqr * n^(-1/3))
    #
    bin_width = np.round(2 * iqr(counts_arr) / np.power(len(counts_arr), 1./3.))
    # calculate number of bins
    num_bins = np.round(int((np.max(counts_arr) - np.min(counts_arr)) / bin_width))

    # get histogram
    hist_counts, hist_bins = np.histogram(counts_arr, int(num_bins * 2.5))

    # find peaks
    # set minimum distance between peaks as half typical ion counts (convert to bins)
    _peak_dist = np.round((90./1.5) / bin_width)
    peaks, props = find_peaks(hist_counts, distance=_peak_dist)

    # assume first and second peaks are the counts for background and 1 ion fluorescence, respectively
    # todo: handle exception case where fewer than 2 peaks
    counts_bgr, counts_signal = hist_bins[peaks[: 2]]

    # use roos thesis minimum error thresholding to calculate appropriate thresholds
    counts_threshold = counts_signal / np.log(1 + counts_signal/counts_bgr)

    # guess number of ions as num_peaks - 1
    num_ions = len(peaks) - 1

    return counts_signal, counts_bgr, counts_threshold, num_ions


def groupBy(dataset, column_num=0, reduce_func=lambda x:x):
    """
    Groups a 2-D array by a given column.
    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # todo: add argument which is a function that can act on grouped result
    # ensure dataset is a numpy array for ease of use
    dataset = np.array(dataset)

    # sort array by given column number (necessary for itertools.groupby)
    dataset = dataset[np.argsort(dataset[:, column_num])]

    # group same values into dict
    dataset_processed = {
        key: reduce_func(np.array([np.delete(val, column_num, 0) for val in group]))
        for key, group in groupby(dataset, lambda arr: arr[column_num])
    }

    return dataset_processed


def importDatasetArtiq(dataset, *args, **kwargs):
    pass
