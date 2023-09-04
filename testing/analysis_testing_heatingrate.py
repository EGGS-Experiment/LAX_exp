"""
Analysis Testing - HeatingRate

Test fitting for Heating Rate
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

from LAX_exp.experiments.SidebandCooling import SidebandCooling
from LAX_exp.experiments.HeatingRate import HeatingRate

# datafile parameters
directory_path =            '/Users/claytonho/Documents/Research/Data & Analysis/Heating Rate/Datasets'
datafile_key =              '31232'


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
    class ExperimentWrapper(SidebandCooling):

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


    # process results for user
    rsb_params = storage_dict['fit_params_rsb']
    bsb_params = storage_dict['fit_params_bsb']
    time_fit_s = parameters['time_readout_pipulse_us']
    print(rsb_params)
    print(bsb_params)

    # separate RSB and BSB
    def split(arr, cond):
        return [arr[cond], arr[~cond]]
    results_rsb, results_bsb = split(exp_res, exp_res[:, 0] < np.mean(exp_res[:, 0]))


    # sinc fit function
    def fit_func(x, a, b, c):
        """
        todo: document arguments
        """
        return ((a**2. / (a**2. + (x - b)**2.)) * np.sin((np.pi * time_fit_s) * (a**2. + (x - b)**2.)**0.5)**2. + c)

    # plot RSB
    rsb_fit = fit_func(results_rsb[:, 0], *rsb_params)
    plt.plot(*(results_rsb.transpose()), 'x')
    plt.plot(*(results_rsb.transpose()))
    plt.plot(results_rsb[:, 0], rsb_fit)
    plt.show()

    # plot BSB
    bsb_fit = fit_func(results_bsb[:, 0], *bsb_params)
    plt.plot(*(results_bsb.transpose()), 'x')
    plt.plot(*(results_bsb.transpose()))
    plt.plot(results_bsb[:, 0], bsb_fit)
    plt.show()


except Exception as e:
    print("Error during testing: {}".format(e))
