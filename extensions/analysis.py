"""
LAX.extensions.analysis

Contains helpful/commonly used analysis modules.
"""

# __all__ = []


# necessary imports
import numpy as np


'''
Dataset Processing
'''
def groupBy(dataset, combine=False):
    """
    Groups a 2-D array by its first column.


    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # ensure dataset is a numpy array for ease of use
    dataset = np.array(dataset)

    # get sorted list of unique x-values (i.e. first column)
    x_values = sorted(set(dataset[:, 0]))

    # group results by x-values
    dataset_processed = {
        x_val: []
        for x_val in x_values
    }
    for x_val, y_values in dataset:
        dataset_processed[x_val].append(y_values)

    # combine y_values into a single list
    if combine:
        for key, value in dataset_processed.items():
            dataset_processed[key] = sum(value, [])

    return dataset_processed

def importDatasetArtiq(dataset, *args, **kwargs):
    pass


'''
Photons
'''
def calculateThreshold(signal, noise, num_ions=1):
    """
    Calculate the discrimination level.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # todo: account for multiple ions
    return signal / np.log(1 + signal / noise)

def discriminateCounts(counts, threshold):
    """
    Process a list of counts into a probability value.
    # todo: accounts for multiple ions

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    # ensure threshold argument is a list
    if not isinstance(threshold, (list, np.ndarray)):
        threshold = [threshold]
    else:
        raise Exception("Error: invalid threshold type.")

    # create holding array
    res = np.zeros(len(counts))

    # discriminate counts
    for cutoff in threshold:
        res += np.where(counts > cutoff, 1, 0)

    return res


'''
Fitting
'''
def fitRabiFlopping(data):
    """
    Fit Rabi Flopping population data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass

def fitSidebandCooling(data):
    """
    Fit sideband cooling spectrum data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass

def fitTemperature(data):
    """
    Fit temperature measurement data.

    Arguments:
        ***todo

    Returns:
        ***todo
    """
    pass
