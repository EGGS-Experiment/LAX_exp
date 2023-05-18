"""
LAX.analysis.processing

Contains helpful/commonly used modules for processing datasets.
"""

# __all__ = []


# necessary imports
import numpy as np


'''
Dataset Processing
'''
def groupBy(dataset, column=1, combine=False):
    """
    Groups a 2-D array by a given column.


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
