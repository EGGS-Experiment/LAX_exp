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
directory_path =            '/Users/claytonho/Documents/Research/Data & Analysis/Parametric/Datasets/2023_12_testing'
datafile_key =              '43684'

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
    # print("Desired file found:\t{}".format(directory_path))


    # extract results and parameters from experiment dataset
    with h5py.File(desired_filepath, 'r') as file:
        results = np.array(file['datasets']['results']).copy()
        parameters = dict(file['arguments'].attrs).copy()
        system = dict(file['system'].attrs).copy()

    # print(list(parameters.keys()))
    # print(system.keys())


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
    # for i in range(2, 5):
    #     im = plt.imshow(th1[:, :, i], interpolation='bicubic')
    #     cbar = plt.colorbar(im)
    #     plt.show()


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
            fitmintmp(ind) for ind in range(len(th1[0, :, 0]))
    }

    # save and print results
    # for key, val in yz0.items():
    #     print('\tFreq. (kHz): {:.2f}\tOpt. Voltage (V): {:.2f}'.format(key, val))

    # print summary statistics
    yz1 = list(yz0.values())
    # print("\nMean (V):\t\t{:.3f}\nMedian (V):\t\t{:.3f}\nStd (V):\t\t{:.3f}".format(
    #     np.mean(yz1),
    #     np.median(yz1),
    #     np.std(yz1))
    # )

    from scipy.optimize import lsq_linear
    def complexLinearFit2(dataset):
        vectorY = dataset[:, 1]
        matrixA = np.array([np.ones(len(vectorY)), dataset[:, 0]]).transpose()
        res = lsq_linear(matrixA, vectorY)
        b_fit_re, b_fit_im = (res.x[0].real, res.x[0].imag)
        m_fit_re, m_fit_im = (res.x[1].real, res.x[1].imag)
        voltage_optimal = - (b_fit_re * m_fit_re + b_fit_im * m_fit_im) / (m_fit_re ** 2. + m_fit_im ** 2.)
        min_pos = res.x[0] + voltage_optimal * res.x[1]

        # tmp remove
        aa_err = np.array([np.sum(np.abs(resids)) / np.sqrt(len(resids)) for resids in aa1])
        # tmp remove
        return np.array([*res.x, min_pos, voltage_optimal, aa_err])
    def fitmintmp2(col_num):
        _restmp = th1[:, col_num, [1, 2, 3]]
        _restmp = np.array([
            _restmp[:, 0],
            _restmp[:, 1] * np.exp(1.j * _restmp[:, 2])
        ], dtype='complex128').transpose()
        return complexLinearFit2(_restmp)

    kk0 = {
        th1[0, ind, 0]*1000.: fitmintmp2(ind)
        for ind in range(len(th1[0, :, 0]))
    }
    kk1 = np.array(list(kk0.values()))
    kk20 = np.array([[np.abs(val), np.angle(val) / np.pi] for val in kk1[:, 0]])
    kk21 = np.array([[np.abs(val), np.angle(val) / np.pi] for val in kk1[:, 1]])
    kk22 = np.array([[np.real(val), np.imag(val) / np.pi] for val in kk1[:, 2]])

    aa0 = kk1[:, 3]
    aa1 = kk1[:, 4]



    # # plot extracted minima (in complex plane)
    # for val in kk22:
    #     plt.plot(*val, 'o')
    #     print(val)
    # plt.plot(*kk22.transpose())
    # plt.show()

    # plt.plot(th1[0, :, 0], np.abs(kk1[:,2]))
    plt.show()

    # plot voltage sweeps in complex plane for each frequency
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.grid(True, which='both')
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    for ind in range(len(th0[0, :, 0])):
        if ind == 7:
            scde = np.array([[np.real(val[0] * np.exp(1.j * val[1])), np.imag(val[0] * np.exp(1.j * val[1]))]
                             for val in th1[:, ind, [2, 3]]])
            # print(scde)
            # ax.plot(*scde.transpose(), marker='x', alpha=0.3)
    # plt.show()

    fixed_channel = 'H Shim' if parameters['dc_micromotion_channel'] == 'V Shim' else 'V Shim'
    print('RID: {:s}'.format(datafile_key))
    print('\tParameters:\n\t\tMode (kHz):\t\t{:.2f}\n\t\tFixed Channel:\t{:s}\n\t\tVoltage:\t\t{:.2f}'.format(
        np.mean(th1[:, 0, 0]) * 1000,
        fixed_channel,
        system['dc_voltage_{:s}_v'.format(fixed_channel)])
    )
    print("\tY-Intercept:\n\t\tMean:\t\t{}\n\t\tMedian:\t\t{}\n\t\tStd:\t\t{}".format(
        np.mean(kk20, axis=0),
        np.median(kk20, axis=0),
        np.std(kk20, axis=0))
    )
    print("\tSlope:\n\t\tMean:\t\t{}\n\t\tMedian:\t\t{}\n\t\tStd:\t\t{}".format(
        np.mean(kk21, axis=0),
        np.median(kk21, axis=0),
        np.std(kk21, axis=0))
    )
    print("\tMin Ampl:\n\t\tMean (%):\t\t{:.3f}\n\t\tMedian (%):\t\t{:.3f}\n\t\tStd (%):\t\t{:.3f}".format(
        np.mean(np.abs(kk1[:, 2]))*100.,
        np.median(np.abs(kk1[:, 2]))*100.,
        np.std(np.abs(kk1[:, 2]))*100.)
    )
    print("\tVoltage ({:s}):\n\t\tMean (V):\t\t{:.3f}\n\t\tMedian (V):\t\t{:.3f}\n\t\tStd (V):\t\t{:.3f}".format(
        parameters['dc_micromotion_channel'],
        np.mean(yz1),
        np.median(yz1),
        np.std(yz1))
    )

except Exception as e:
    print("Error during testing: {}".format(repr(e)))
