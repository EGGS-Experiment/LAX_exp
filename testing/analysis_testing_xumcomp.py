"""
Analysis Testing - Micromotion Compensation

Test micromotion compensation.
"""
'''
IMPORTS
'''
# general imports
import os
import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt

# ARTIQ & LAX imports
from artiq.experiment import *
from LAX_exp.analysis import *
from LAX_exp.extensions import *

# from LAX_exp.experiments.parametric.ParametricSweep import ParametricSweep

# datafile parameters
directory_path =            '/Users/claytonho/Documents/Research/Data & Analysis/Parametric/Datasets'
datafile_key =              '43722'


'''
MAIN SEQUENCE
'''
try:
    # search given directory path for file with desired key
    def find_file(search_directory, search_key):
        for root, dirs, files in os.walk(search_directory):
            for filename in files:
                if search_key in filename:
                    return os.path.join(root, filename)
        raise Exception("Error: filekey {} not found in {}".format(search_directory, search_key))

    desired_filepath = find_file(directory_path, datafile_key)
    print("Desired file found:\t{}".format(directory_path))


    # extract results and parameters from experiment dataset
    with h5py.File(desired_filepath, 'r') as file:
        results = np.array(file['datasets']['results']).copy()
        parameters = dict(file['arguments'].attrs).copy()


    from itertools import groupby
    def groupBy2(dataset, column_num=0, reduce_func=lambda x: x):
        # ensure dataset is a numpy array for ease of use
        dataset = np.array(dataset)

        # sort array by given column number (necessary for itertools.groupby)
        dataset = dataset[np.argsort(dataset[:, column_num])]

        # group same values into dict
        dataset_processed = [
            [key, reduce_func(np.array([val for val in group]))]
            for key, group in groupby(dataset, lambda arr: arr[column_num])
        ]

        return dataset_processed

    # group in 2D
    th0 = np.array([tmp[1] for tmp in groupBy2(results, column_num=1)])[::-1]
    th1 = np.array([x[np.argsort(x[:, 0])] for x in th0])

    # plot correlated amplitude, correlated phase, and count rate
    for i in range(2, 5):
        im = plt.imshow(th1[:, :, i], interpolation='bicubic')
        cbar = plt.colorbar(im)
        plt.show()


    '''
    MAGICAL FITTING
    '''
    def fitmintmp(col_num):
        _restmp = th1[:, col_num, [1, 2, 3]]
        _restmp = np.array([
            _restmp[:, 0],
            _restmp[:, 1] * np.exp(1.j * _restmp[:, 2])
        ], dtype='complex128').transpose()
        return complexLinearFitMinimize(_restmp)

    yz0 = {
        th1[0, ind, 0]*1000.:
            fitmintmp(ind) for ind in range(len(th1))
    }
    yz1 = list(yz0.values())
    print(yz1)
    print("\nMean (V):\t\t{:.3f}\nMedian (V):\t\t{:.3f}\nStd (V):\t\t{:.3f}".format(
        np.mean(yz1),
        np.median(yz1),
        np.std(yz1))
    )


except Exception as e:
    print("Error during testing: {}".format(e))
