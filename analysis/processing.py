"""
LAX.analysis.processing

Contains helpful/commonly used modules for processing datasets.
"""

__all__ = ['findThresholdScikit', 'findThresholdPeaks', 'groupBy', 'processDataset2D']


# necessary imports
import numpy as np
from itertools import groupby

# from scipy.stats import iqr
from scipy.signal import find_peaks
# from skimage.filters import threshold_otsu, threshold_multiotsu, threshold_minimum, threshold_yen, threshold_isodata, threshold_triangle
from skimage.filters import threshold_multiotsu, threshold_minimum


'''
Thresholding
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
        bin_width = 3.49 * np.std(counts_arr) / (len(counts_arr) ** (1. / 3.))
        # calculate number of bins
        num_bins = round((np.max(counts_arr) - np.min(counts_arr)) / bin_width)

    # guess number of ions (i.e. thresholding classes) if no value is provided
    if num_ions is None:
        # histogram counts
        hist_counts, hist_bins = np.histogram(counts_arr, int(num_bins * 1.8))
        # guess number of ions using heuristic values (~20 background counts, ~150 max counts per ion)
        num_ions = round((np.max(hist_bins) - 20.) / 150.)

    # tmp remove
    if num_ions == 0: num_ions = 1
    # tmp remove

    # start with minimum error thresholding since it's most likely to recognize the signal/background threshold
    thresh_values_start =  np.array([threshold_minimum(counts_arr)])
    # use multi-otsu thresholding to get threshold values in case of num_ions > 1
    thresh_multiotsu_values = threshold_multiotsu(counts_arr, classes=num_ions+1, nbins=num_bins)

    # ensure duplicate thresholds are not added to list
    for thresh_val in thresh_multiotsu_values:
        # only recognize value if different from existing thresholds by thresh_dist
        if np.all(np.abs(thresh_values_start - thresh_val) > thresh_dist):
            thresh_values_start = np.append(thresh_values_start, thresh_val)


    # # create list of threshold functions
    # threshold_functions = [threshold_isodata, threshold_yen, threshold_triangle]
    #
    # # ensure duplicate thresholds are not added to list
    # for thresh_func in threshold_functions:
    #     # get threshold value
    #     thresh_val = thresh_func(counts_arr, nbins=num_bins)
    #
    #     # only recognize value if different from existing thresholds by thresh_dist
    #     if np.all((thresh_values - thresh_val) > thresh_dist):
    #         thresh_values = np.append(thresh_values, thresh_val)

    # tmp remove
    print('\t\t\tthresholds: {}'.format(np.sort(thresh_values_start)[: num_ions]))
    # tmp remove

    # return sorted threshold values
    return np.sort(thresh_values_start)[: num_ions]


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
    bin_width = 3.49 * np.std(counts_arr) / (len(counts_arr) ** (1. / 3.))
    # calculate number of bins
    num_bins = round((np.max(counts_arr) - np.min(counts_arr)) / bin_width)

    # get histogram
    hist_counts, hist_bins = np.histogram(counts_arr, int(num_bins * 2.5))
    # plt.stairs(hist_counts, hist_bins)

    # find peaks
    # set minimum distance between peaks as half typical ion counts (convert to bins)
    _peak_dist = round((90./1.5) / bin_width)
    peaks, props = find_peaks(hist_counts, distance=_peak_dist)

    # assume first and second peaks are the counts for background and 1 ion fluorescence, respectively
    # todo: handle exception case where fewer than 2 peaks
    counts_bgr, counts_signal = hist_bins[peaks[: 2]]

    # use roos thesis minimum error thresholding to calculate appropriate thresholds
    counts_threshold = counts_signal / np.log(1 + counts_signal/counts_bgr)

    # guess number of ions as num_peaks - 1
    num_ions = len(peaks) - 1

    return counts_signal, counts_bgr, counts_threshold, num_ions


'''
Dataset Processing
'''
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


def processDataset2D(dataset):
    """
    todo: document

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # create data structures for processing
    probability_vals =      np.zeros(len(dataset))

    # calculate fluorescence detection threshold
    threshold_list =        findThresholdScikit(dataset[:, 1])
    for threshold_val in threshold_list:
        probability_vals[np.where(dataset[:, 1] > threshold_val)] += 1.
    # normalize probabilities and convert from D-state probability to S-state probability
    dataset[:, 1] =         1. - probability_vals / len(threshold_list)

    # process dataset into x, y, with y being averaged probability
    dataset_processed =     groupBy(dataset, column_num=0, reduce_func=np.mean)
    dataset_processed =     np.array([list(dataset_processed.keys()),
                                      list(dataset_processed.values())]).transpose()
    return dataset_processed
