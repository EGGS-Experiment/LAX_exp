"""
LAX.analysis.processing

Contains helpful/commonly used modules for processing datasets.
"""

__all__ = ['findThresholdScikit', 'findThresholdPeaks',
           'groupBy', 'groupBy2',
           'processFluorescence2D', 'extract_ratios', 'extract_sidebands_freqs', 'convert_ratios_to_coherent_phonons',
           'convert_ratios_to_squeezed_phonons', 'process_laser_scan_results']


# necessary imports
import numpy as np
from itertools import groupby

# from scipy.stats import iqr
from scipy.signal import find_peaks
# from skimage.filters import threshold_otsu, threshold_multiotsu, threshold_minimum, threshold_yen, threshold_isodata, threshold_triangle
from skimage.filters import threshold_multiotsu, threshold_minimum
from scipy.special import factorial
from scipy.interpolate import interp1d

from LAX_exp.extensions.conversions import *


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
    thresh_values_start = np.array([threshold_minimum(counts_arr)])
    # use multi-otsu thresholding to get threshold values in case of num_ions > 1
    thresh_multiotsu_values = threshold_multiotsu(counts_arr, classes=num_ions + 1, nbins=num_bins)

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

    # print result to log
    # print('\tThresholds: {}'.format(np.sort(thresh_values_start)[: num_ions]))
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
    _peak_dist = round((90. / 1.5) / bin_width)
    peaks, props = find_peaks(hist_counts, distance=_peak_dist)

    # assume first and second peaks are the counts for background and 1 ion fluorescence, respectively
    # todo: handle exception case where fewer than 2 peaks
    counts_bgr, counts_signal = hist_bins[peaks[: 2]]

    # use roos thesis minimum error thresholding to calculate appropriate thresholds
    counts_threshold = counts_signal / np.log(1 + counts_signal / counts_bgr)

    # guess number of ions as num_peaks - 1
    num_ions = len(peaks) - 1

    return counts_signal, counts_bgr, counts_threshold, num_ions


'''
Dataset Processing
'''

def groupBy(dataset, column_num=0, reduce_func=lambda x: x):
    """
    Groups a 2-D array by a given column.
    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # ensure dataset is a numpy array for ease of use
    dataset = np.array(dataset)

    # sort array by given column number (necessary for itertools.groupby)
    dataset = dataset[np.argsort(dataset[:, column_num])]

    # group same values into dict
    dataset_processed = {
        # todo: array-ize the np.delete step
        key: reduce_func(np.array([np.delete(val, column_num, 0) for val in group]))
        for key, group in groupby(dataset, lambda arr: arr[column_num])
    }

    return dataset_processed

def groupBy2(dataset, column_nums=0, reduce_func=lambda x: x):
    """
    Groups a 2-D array by a given column.
    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # ensure dataset is a numpy array for ease of use
    dataset = np.array(dataset)

    # ensure column_num is a list
    if type(column_nums) not in (list, np.ndarray):
        column_num = [column_nums]

    def _group(_dataset, _col_num):
        # sort array by given column number (necessary for itertools.groupby)
        _dataset = _dataset[np.argsort(_dataset[:, _col_num])]

        # group same values into list
        return [
            [reduce_func(np.array([val for val in subgroup]))]
            for key, subgroup in groupby(_dataset, lambda arr: arr[_col_num])
        ]

    # for index in column_nums:
    # sort array by given column number (necessary for itertools.groupby)
    dataset = dataset[np.argsort(dataset[:, column_num])]

    # group same values into dict
    dataset_processed = [
        [key, reduce_func(np.array([val for val in group]))]
        for key, group in groupby(dataset, lambda arr: arr[column_num])
    ]

    return dataset_processed

def processFluorescence2D(dataset):
    """
    todo: document

    Args:
        ***todo

    Returns:
        ***todo
    """
    # create data structures for processing
    probability_vals = np.zeros(len(dataset))

    # calculate fluorescence detection threshold
    threshold_list = findThresholdScikit(dataset[:, 1])
    for threshold_val in threshold_list:
        probability_vals[np.where(dataset[:, 1] > threshold_val)] += 1.
    # normalize probabilities and convert from D-state probability to S-state probability
    dataset[:, 1] = 1. - probability_vals / len(threshold_list)

    # process dataset into x, y, with y being averaged probability
    dataset_processed = groupBy(dataset, column_num=0, reduce_func=np.mean)
    dataset_processed = np.array([list(dataset_processed.keys()),
                                  list(dataset_processed.values())]).transpose()
    return dataset_processed


"""
EGGS HEATING FUNCTIONALITY
"""

def extract_ratios(dataset: np.array,
                   sorting_col_num: int, counts_col_num: int, readout_col_num: int,
                   reps: int, sub_reps: int):
    """
    Calculate the rsb/bsb ratios of a dataset.

    Arguments:
        dataset: dataset to be analyzed
        sorting_col_num: column of the dataset contain the frequencies that were scanned (sideband, carrier, etc.)
        counts_col_num: column number of the dataset containing fluorescence counts
        readout_col_num: column number of the dataset countaining readout frequencies
        reps: number of repetitions performed at each experimental point
        sub_reps: number of sub-reps performed for each repetition

    Returns:
        ratios: rsb/bsb ratios
        probs_rsb: the rsb excitation probability
        probs_bsb: the rsb excitation probability
        std_rsb: standard deviation for the rsb excitation probability
        std_bsb: standard deviation for the bsb excitation probability
        scanning_freqs_MHz_unique: frequencies we scan over
    """
    dataset_sorted = dataset[np.argsort(dataset[:, sorting_col_num]), :]
    scanning_freqs = dataset_sorted[:, sorting_col_num]
    scanning_freqs_unique = np.unique(scanning_freqs)
    readout_freqs_sorted = np.array(dataset_sorted[:, readout_col_num])
    counts = np.array(dataset_sorted[:, counts_col_num])


    if np.array_equal(scanning_freqs, readout_freqs_sorted):
        scanning_freqs_MHz_unique = scanning_freqs_unique * (2 * 2.32830644e-7)
    else:
        scanning_freqs_MHz_unique = scanning_freqs_unique * 1e-6

    readout_freqs_MHz_sorted = readout_freqs_sorted* (2 * 2.32830644e-7)
    probs = np.zeros(len(counts))
    guess_Ca_carrier_MHz = np.mean(np.unique(readout_freqs_MHz_sorted))

    # determine thresholds
    threshold_list = findThresholdScikit(counts)
    for threshold_val in threshold_list:
        probs[np.where(counts > threshold_val)] += 1.

    normalized_probs = 1. - probs / len(threshold_list)
    normalized_probs_rsb = normalized_probs[guess_Ca_carrier_MHz > readout_freqs_MHz_sorted]
    probs_rsb = np.mean(normalized_probs_rsb.reshape(-1, sub_reps * reps), 1)

    std_rsb = np.std(np.reshape(normalized_probs_rsb, (-1, reps * sub_reps)), 1)  / np.sqrt(reps * sub_reps)

    normalized_probs_bsb = normalized_probs[guess_Ca_carrier_MHz < readout_freqs_MHz_sorted]
    probs_bsb = np.mean(np.reshape(normalized_probs_bsb, (-1, reps * sub_reps)), 1)
    std_bsb = np.std(np.reshape(normalized_probs_bsb, (-1, reps * sub_reps)), 1) / np.sqrt(reps * sub_reps)

    ratios = np.divide(probs_rsb, probs_bsb + 1e-7)
    return ratios, probs_rsb, probs_bsb, std_rsb, std_bsb, scanning_freqs_MHz_unique

def extract_sidebands_freqs(readout_freqs_MHz):
    """
    Split readout frequencies into red and blue sidebands - should already be sorted and unique

    Args:
        readout_freqs_MHz: sorted and unique list of frequencies to readout for sideband analysis

    Returns:
        rsb_freqs: frequencies for readout of rsb
        bsb_freqs: frequencies for readout of bsb
        guess_Ca_carrier_MHz: expected carrier frequency
    """
    guess_Ca_carrier_MHz = np.mean(readout_freqs_MHz)
    rsb_freqs = readout_freqs_MHz[guess_Ca_carrier_MHz > readout_freqs_MHz]
    bsb_freqs = readout_freqs_MHz[guess_Ca_carrier_MHz < readout_freqs_MHz]
    return rsb_freqs, bsb_freqs, guess_Ca_carrier_MHz


"""
Functions for Coherent States
"""

def convert_ratios_to_coherent_phonons(ratios: np.array) -> np.array:
    """
    Convert rsb/bsb ratios to number of phonons for a coherent state

    Argus:
        ratios: rsb/bsb ratios from sidebands

    Returns:
        phonons: phonon count of coherent state
    """
    ratios[ratios < 0] = 0
    ratios[ratios > .8] = .8

    nbars = np.linspace(0, 2, 2001)
    coherent_ratios = np.zeros(len(nbars))
    for idx, nbar in enumerate(nbars):
        coherent_ratios[idx] = prob_rsb_coherent(nbar) / prob_bsb_coherent(nbar)

    interp_func = interp1d(coherent_ratios, nbars)
    phonons = interp_func(ratios)
    return phonons

def coherent_state_amp(nbar, n):
    """
    Determine probability amplitudes of a coherent state
    Arguemnts:
        nbar: average phonon number
        n: phonon number
    Returns: coherent state amplitude


    todo: support taking array of nbar (currently only supports n as an array but nbar needs to be an int/float)
    """
    return np.multiply(np.exp(-np.abs(nbar) / 2), np.power(np.sqrt(nbar), n) / np.sqrt(factorial(n)))

def prob_bsb_coherent(nbar):
    """
    Determine blue sideband excitation probability for a coherent state

    Args:
        nbar: average phonon number

    Returns:
        blue sideband excitation prbability
    """
    n = np.arange(0, 100)
    return 1 - 1 / 2 * np.sum((1 + np.cos(np.pi * np.sqrt(n + 1))) * np.abs(coherent_state_amp(nbar, n)) ** 2)

def prob_rsb_coherent(nbar):
    """
    Determine red sideband excitation probability for a coherent state

    Args:
        nbar: average phonon number

    Returns:
        red sideband excitation probability
    """
    n = np.arange(1, 100)
    return 1 - np.abs(coherent_state_amp(nbar, 0)) ** 2 - 1 / 2 * np.sum(
        (1 + np.cos(np.pi * np.sqrt(n))) * np.abs(coherent_state_amp(nbar, n)) ** 2)


"""
Functions for Squeezed States
"""
def convert_ratios_to_squeezed_phonons(ratios: np.array) -> np.array:
    """
    Convert rsb/bsb ratios to number of phonons for a squeeze state

    Arguments:
        ratios: rsb/bsb ratios from sidebands

    Returns:
        phonons: phonon count of squeezed state
    """
    ratios[ratios < 0] = 0
    ratios[ratios > .8] = .8

    rs = np.linspace(0, 2.0, 2001)
    squeeze_ratios = np.zeros(len(rs))
    for idx, r in enumerate(rs):
        squeeze_ratios[idx] = prob_rsb_squeeze(r) / prob_bsb_squeeze(r)

    interp_func = interp1d(squeeze_ratios, np.sinh(rs)**2)
    phonons = interp_func(ratios)
    return phonons

def squeeze_state_population(r,n):
    if isinstance(n, int):
        if n<10:
            return (np.tanh(r)**(2*n))/(np.cosh(r))*(factorial(2*n))/((2**n) * factorial(n))**2
        else:
            return (np.tanh(r)**(2*n))/(np.cosh(r))*1/np.sqrt(np.pi*n)
    low_n = n[n<10]
    low =  (np.tanh(r)**(2*low_n))/(np.cosh(r))*(factorial(2*low_n))/((2**low_n) * factorial(low_n))**2
    high_n = n[n>=10]
    high =  (np.tanh(r)**(2*high_n))/(np.cosh(r))*1/np.sqrt(np.pi*high_n)
    return np.concatenate((low,high))

def prob_rsb_squeeze(r):
    n = 2*np.arange(1, 35)
    return 1 - squeeze_state_population(r,0) - 1 / 2 * np.sum((1 + np.cos(np.pi * np.sqrt(n))) * squeeze_state_population(r,n/2))

def prob_bsb_squeeze(r):
    n = 2*np.arange(0, 35)
    return 1 - 1 / 2 * np.sum((1+np.cos(np.pi * np.sqrt(n + 1))) * squeeze_state_population(r,n/2))


"""
Laser Scan Functionality
"""
def process_laser_scan_results(results, time_us):
    # todo: move to use processFluorescence2D
    # create data structures for processing
    results_tmp =           np.array(results)[:, :2]
    probability_vals =      np.zeros(len(results_tmp))
    counts_arr =            np.array(results_tmp[:, 1])

    # convert x-axis (frequency) from frequency tuning word (FTW) to MHz
    results_tmp[:, 0] *=    1.e3 / 0xFFFFFFFF

    # calculate fluorescence detection threshold
    threshold_list =        findThresholdScikit(results_tmp[:, 1])
    for threshold_val in threshold_list:
        probability_vals[np.where(counts_arr > threshold_val)] += 1.
    # normalize probabilities and convert from D-state probability to S-state probability
    results_tmp[:, 1] =     1. - probability_vals / len(threshold_list)

    # process dataset into x, y, with y being averaged probability
    results_tmp =           groupBy(results_tmp, column_num=0, reduce_func=np.mean)
    results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()


    try:
        # calculate peak criteria from data
        # todo: somehow relate peak height to shot noise (i.e. 1/sqrt(N))
        # todo: maybe set min peak width of at least 2 points (? not sure if good idea)
        # _peak_height =          np.power(repetitions, -0.5)
        _peak_height =          0.2
        _peak_thresh =          0.05
        # peak distance criteria is set as ~8 kHz between points
        _peak_dist =            int(4e-3 / (results_tmp[1, 0] - results_tmp[0, 0]))

        # calculate peaks from data and extract values
        from scipy.signal import find_peaks
        peaks, props =          find_peaks(results_tmp[:, 1], height=_peak_height, distance=_peak_dist)
        peak_vals =             results_tmp[peaks]

        # fit sinc profile to results (only in the case of one peak)
        if len(peaks) == 1:
            # get index step size in frequency (mhz)
            step_size_mhz = np.mean(results_tmp[1:, 0] - results_tmp[:-1, 0])
            freq_sinc_mhz = 1. / time_us

            # get points +/- 6x the 1/f time for sinc fitting
            num_points_sinc = round(6. * freq_sinc_mhz / step_size_mhz)
            index_peak_center = peaks[0]
            index_min = max(0, index_peak_center - num_points_sinc)
            index_max = min(index_peak_center + num_points_sinc, len(results_tmp))
            points_tmp = results_tmp[index_min: index_max]

            # fit sinc profile and replace spectrum peak with fitted value
            # note: division by 2 accounts for conversion between AOM freq. and abs. freq.
            from LAX_exp.analysis.fitting import fitSinc
            fitter = fitSinc()
            fit_sinc_params, _ = fitter.fit(points_tmp, time_us / 2.)
            peak_vals[0, 0] = fit_sinc_params[1]

    except Exception as e:
        peak_vals = []

    return peak_vals, results_tmp
