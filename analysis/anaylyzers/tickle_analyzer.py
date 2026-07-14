# load from external libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# local imports
from base.dataloader import DataLoader
from base.artiq_conversions import *
from base.units import *
from base.processing.processing import convert_and_sort
from base.processing.processing import process_multi_ion_state_vals
from LAX_exp.analysis.processing import findThresholdScikit
from base.fitting.cat_fitting_functions import *
from base.fitting.fitting_functions import *
from base.fitting.motional_fitting import *
from base.conversions.physics_conversions import *

from numpy import array, int32, int64, zeros
from numpy import shape as array_shape


class TickleAnalyzer:

    def __init__(self, dataloader, rid, detuning_idx):

        # grab the data
        self.dataloader = dataloader
        self.rid = rid
        self.reps = dataloader.arguments[rid]['repetitions']
        self.detuning_idx = detuning_idx

    def analyze_tickle_experiment(self):
        results, num_states = self._process_results()

        # preparation for conversion to phonons
        alphas = linspace(0, 4, 100)
        nbars = [alpha_to_nbar(alpha) for alpha in alphas]

        phonon_conversions = [phonon_conversion(alpha) for alpha in alphas]
        time_tickle_us = self.dataloader.arguments[self.rid]['time_heating_us']

        detunings = results['detuning']
        population_vals = results['populations']
        population_errs = results['populations_err']

        # phonons = interp(population_vals[:,0], phonon_conversions, nbars)
        #
        # # from scipy.differentiate import derivative
        # # phonon_errs = derivative(p, x_point)
        # phonon_errs = [np.nan]*len(detunings)
        #
        fit_x = np.linspace(np.min(detunings), np.max(detunings), 1000)

        p = population_vals[:, 0]
        p_err = population_errs[:, 0]

        # np.interp needs increasing x-values
        x = np.asarray(phonon_conversions)
        y = np.asarray(nbars)

        order = np.argsort(x)
        x = x[order]
        y = y[order]

        phonons = np.interp(p, x, y)

        # local derivative dnbar / dpopulation
        slopes = np.gradient(y, x)
        local_slopes = np.interp(p, x, slopes)

        phonon_errs = np.abs(local_slopes) * p_err

        '''guess parameters'''
        max_phonon = max(phonons)
        peak_freq = detunings[np.argmax(phonons)]
        linewidth_scalar = time_tickle_us/1e3
        contrast_loss = 0
        base_p0 = [max_phonon, peak_freq, linewidth_scalar, contrast_loss]


        bounds = (
            [0.0, -np.inf, 0.0, 0.0, 0.0],
            [10, np.inf, np.inf, 1],
        )

        try:
            popt, pcov = curve_fit(sinc_fit, detunings, phonons,
                                       p0=base_p0, maxfev=15000)

            perr = np.sqrt(np.diag(pcov))

            fit_y = sinc_fit(fit_x, *popt)

        except:
            print(f'Could not find fit')
            fit_y= [np.nan] * len(fit_x)

        raw_data = {'x': detunings,
                    'y': phonons,
                    'yerr': phonon_errs,
                    'legend_labels': 'Phonons'
        }

        fittting_results = {
            'fit_x': fit_x,
            'fit_y': fit_y,
            'popt': popt,
            'pcov': pcov,
        }

        return raw_data, fittting_results, num_states


    def _process_results(self):
        # get results
        results_tmp = self.dataloader.results[self.rid]

        counts_arr = array(results_tmp[:, 1])
        detuning_arr = ftw_to_frequency_khz(array(results_tmp[:, self.detuning_idx]))

        detuning_vals = np.unique(detuning_arr)

        # calculate fluorescence detection threshold
        threshold_list = np.sort(findThresholdScikit(counts_arr))
        num_states = len(threshold_list) + 1

        results_storer = {}

        count_states = np.digitize(counts_arr, threshold_list)

        # group all shots by identical detuning
        detunings, detuning_idx = np.unique(detuning_arr, return_inverse=True)
        population_vals = np.zeros((len(np.unique(detunings)), num_states))
        population_err = np.zeros((len(np.unique(detunings)), num_states))

        for det_idx in range(len(detunings)):
            det_mask = det_idx == detuning_idx
            for state in range(num_states):
                counts_masked = count_states[det_mask] == state
                population_vals[det_idx, state] = np.mean(counts_masked)
                population_err[det_idx, state] = np.std(counts_masked) / np.sqrt(len(counts_masked))

        results_storer = {
            "detuning": detunings,
            "populations": population_vals,
            "populations_err": population_err
        }

        return results_storer, num_states