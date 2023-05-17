"""
LAX.analysis.segregation

Contains modules used for segregation analysis of datasets (e.g. thresholding).
"""

# __all__ = []


# necessary imports
import numpy as np



'''
Thresholding/Discrimination
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
    # ensure counts is a numpy array
    counts = np.array(counts)

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

    return np.mean(res)
