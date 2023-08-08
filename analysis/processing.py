"""
LAX.analysis.processing

Contains helpful/commonly used modules for processing datasets.
"""

__all__ = ['findThresholdScikit', 'findThresholdPeaks', 'groupBy']


# necessary imports
import numpy as np

from itertools import groupby

from scipy.stats import iqr
from scipy.signal import find_peaks
from skimage.filters import threshold_otsu, threshold_multiotsu, threshold_minimum, threshold_yen, threshold_isodata, threshold_triangle

# todo: move from np.power to **


'''
Dataset Processing
'''
def findThresholdScikit(counts_arr, thresh_dist=50, num_bins=None, num_ions=None):
    """
    Get the binary discrimination threshold for a dataset
    by using thresholding methods from scikit-image.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # todo: implement error handling somehow

    # calculate num_bins if no value is provided
    if num_bins is None:
        # calculate histogram bin width using freedman diaconis rule (bin width = 2 * iqr * n^(-1/3))
        # bin_width = np.round(2 * iqr(counts_arr) / np.power(len(counts_arr), 1./3.))
        # edit: actually, use scott's normal reference rule instead
        bin_width = np.round(3.49 * np.std(counts_arr) / np.power(len(counts_arr), 1. / 3.))
        # calculate number of bins
        num_bins = np.round(int((np.max(counts_arr) - np.min(counts_arr)) / bin_width))

    # guess number of ions (i.e. thresholding classes) if no value is provided
    if num_ions is None:
        # histogram counts
        hist_counts, hist_bins = np.histogram(counts_arr, int(num_bins * 2.5))
        # get max count bin
        max_counts = np.max(hist_bins)
        # guess number of ions using heuristic values (~20 background counts, ~150 max counts per ion)
        num_ions = int(np.round((max_counts - 20.) / 150.))


    # todo: fix thresholding somehow - want to use minimum error thresholding as the first
    # todo: maybe we can start with minimum threshold value, and if num_ions > 1, then we can start with multiotsu instead
    # use multi-otsu thresholding to get base list of threshold values
    thresh_values = threshold_multiotsu(counts_arr, classes=num_ions+1, nbins=num_bins)
    # create list of threshold functions
    threshold_functions = [threshold_minimum, threshold_isodata, threshold_yen, threshold_triangle]

    # ensure duplicate thresholds are not added to list
    for thresh_func in threshold_functions:
        # get threshold value
        thresh_val = thresh_func(counts_arr, nbins=num_bins)

        # only recognize value if different from existing thresholds by thresh_dist
        if np.all((thresh_values - thresh_val) > thresh_dist):
            thresh_values = np.append(thresh_values, thresh_val)

    # return sorted threshold values
    return np.sort(thresh_values)[: num_ions]


def findThresholdPeaks(counts_arr):
    """
    Get the binary discrimination threshold for a dataset
    by finding histogram peaks and applying minimum error thresholding.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # calculate histogram bin width using freedman diaconis rule (bin width = 2 * iqr * n^(-1/3))
    # bin_width = np.round(2 * iqr(counts_arr) / np.power(len(counts_arr), 1./3.))
    # edit: actually, use scott's normal reference rule instead
    bin_width = np.round(3.49 * np.std(counts_arr) / np.power(len(counts_arr), 1. / 3.))
    # calculate number of bins
    num_bins = np.round(int((np.max(counts_arr) - np.min(counts_arr)) / bin_width))

    # get histogram
    hist_counts, hist_bins = np.histogram(counts_arr, int(num_bins * 2.5))
    # plt.stairs(hist_counts, hist_bins)

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
