"""
Analysis Testing - Laser Scan

Test fitting for LaserScan.
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

from LAX_exp.experiments.diagnostics.LaserScan import LaserScan

# datafile parameters
directory_path =    '/Users/claytonho/Documents/Research/Data & Analysis/Laser Scan/Datasets'
datafile_key =      '37774'


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


    # create experiment object wrapper class to simulate the experiment object
    class ExperimentWrapper(LaserScan):

        def __init__(self, results, parameter_dict, _results_storage={}):
            # instantiate necessary attributes
            self.children = []

            # use results and parameters to instantiate attributes
            self.results = results
            for key, val in parameter_dict.items():
                setattr(self, key, val)

            # redefine set_dataset to squirrel away any data passed to set_dataset
            self.set_dataset = lambda key, val: _results_storage.update({key: val})


    # create experiment instance to wrap around our dataset
    storage_dict = {}
    exp_obj = ExperimentWrapper(results, parameters, storage_dict)
    # call the analyze method of the experiment object to process our results
    exp_res = exp_obj.analyze_experiment()

    plt.plot(*exp_res.transpose())
    plt.show()

    # print(fit_sinc_params)
    # def fit_func(x, a, b, c):
    #     return ((a ** 2. / (a ** 2. + (x - b) ** 2.)) * np.sin(
    #         (np.pi * self.time_qubit_us) * (a ** 2. + (x - b) ** 2.) ** 0.5) ** 2. + c)
    #
    # th1 = fit_func(points_tmp[:, 0], *fit_sinc_params)
    # import matplotlib.pyplot as plt
    # # plt.plot(*results_tmp.transpose())
    # plt.plot(*points_tmp.transpose(), alpha=0.5)
    # plt.plot(points_tmp[:, 0], th1, 'x')
    # plt.show()


except Exception as e:
    print("Error during testing: {}".format(repr(e)))
