# load from external libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# local imports
from LAX_exp.analysis.fitting_functions.cat_fitting_functions import *
from LAX_exp.analysis.artiq_conversions import *
from LAX_exp.analysis.processing import findThresholdScikit

from numpy import array, int32, int64, zeros
from numpy import shape as array_shape
import numpy as np


class CatInterferometerAnalyzer:

    def __init__(self):
        pass

    def analyze_interferometer_experiment(self,
                                          data,
                                          detuning_idx,
                                          time_tickle_us,
                                          enable_ms_gate,
                                          time_cat_bichromatic_us,
                                          time_ms_gate_us,
                                          secular_idx = 0,
                                          use_secular = False):

        data_storer, num_states = self._process_data(data,
                                                     detuning_idx,
                                                     secular_idx = secular_idx,
                                                     use_secular = use_secular)

        try:
            fit_funcs = self._get_fit_funcs(num_states)
        except NotImplementedError as e:
            print("Unable to find fit funcs")
            fit_funcs = None

        if use_secular:
            raw_data = {}
            fitting_results = {}

        for secular_idx in data_storer.keys():
            fit_x_list = []
            fit_y_list = []


            secular = data_storer[secular_idx]['secular']
            detunings = data_storer[secular_idx]['detuning']
            population_vals = data_storer[secular_idx]['populations']
            population_errs = data_storer[secular_idx]['populations_err']

            detuning_list = [detunings for _ in range(num_states)]
            population_list = [population_vals[:, state] for state in range(num_states)]
            population_err_list = [population_errs[:, state] for state in range(num_states)]

            fisher_info = self.get_fisher_info(detuning_list, population_list)
            detuning_list = [detunings for _ in range(num_states +1)]

            fit_x = np.linspace(np.min(detunings), np.max(detunings), 1000)
            fit_x_list = [fit_x for _ in range(num_states)]

            '''guess parameters'''
            d0 = 0.02
            t_guess = time_tickle_us / 1000.0
            alpha0 = 1
            delta_0_guess = 0.0
            phi_guess = 0.0
            # ignores ramp/setup times
            t0_guess = (
                               (time_ms_gate_us if enable_ms_gate else 0.0)
                               + time_cat_bichromatic_us
                       ) / 1000.0

            base_p0 = [d0, alpha0, t_guess, t0_guess, delta_0_guess, phi_guess]

            bounds = (
                [0.0, 0.0, 0.0, 0.0, np.min(detunings), -2 * np.pi],
                [0.5, 20.0, np.inf, np.inf, np.max(detunings), 2 * np.pi],
            )

            if fit_funcs is not None:

                best = self.fit_with_alpha_multistart(
                    detunings,
                    population_vals,
                    fit_funcs,
                    num_states,
                    base_p0,
                    bounds,
                    fit_x,
                )

                fit_y_list = best['fit_y_list']
                popt = best["popt"]
                pcov = best["pcov"]

                print("Best alpha initial guess:", best["alpha0"])
                print("Best fit params:", popt)
                print("RSS:", best["rss"])

            else:
                fit_y_list = [[None] * len(fit_x) for _ in range(num_states)]

            legend_labels = self._make_legend_labels(num_states)

            fisher_info = self.get_fisher_info(fit_x_list, fit_y_list)
            detuning_list = [detunings for _ in range(num_states + 1)]

            ys = np.zeros((num_states + 1, len(detunings)))
            for idx in range(num_states):
                ys[idx] = population_vals[:, idx]
            ys[-1] = [np.nan] * len(detunings)

            errs = np.zeros((num_states + 1, len(detunings)))
            for idx in range(num_states):
                errs[idx] = population_errs[:, idx]
            errs[-1] = [np.nan] * len(detunings)

            fit_xs = np.zeros((num_states + 1, len(fit_x)))
            for idx in range(num_states + 1):
                fit_xs[idx] = fit_x

            fit_ys = np.zeros((num_states + 1, len(fit_x)))
            for idx in range(num_states):
                fit_ys[idx] = fit_y_list[idx]
            fit_ys[-1] = fisher_info

            if use_secular:
                raw_data[secular_idx] = {
                            'x': detuning_list,
                            'y': ys,
                            'yerr': errs,
                            'legend_labels': legend_labels,
                            'secular': secular
                            }
            else:
                raw_data = {
                    'x': detuning_list,
                    'y': ys,
                    'yerr': errs,
                    'legend_labels': legend_labels
                }

            idx_positive_detunings = np.where(fit_x > 0)[0]
            positive_fit_detunings = fit_x[idx_positive_detunings]

            idx_negative_detunings = np.where(fit_x < 0)[0]
            negative_fit_detunings = fit_x[idx_negative_detunings]

            fisher_info_positive_detunings = fisher_info[idx_positive_detunings]
            best_fisher_info_positive_detunings = np.max(fisher_info_positive_detunings)
            best_positive_detuning = positive_fit_detunings[
                np.argmax(fisher_info_positive_detunings)
            ]

            fisher_info_negative_detunings = fisher_info[idx_negative_detunings]
            best_fisher_info_negative_detunings = np.max(fisher_info_negative_detunings)
            best_negative_detuning = negative_fit_detunings[
                np.argmax(fisher_info_negative_detunings)
            ]

            if use_secular:
                fitting_results[secular_idx] = {
                    'fit_x': fit_xs,
                    'fit_y': fit_ys,
                    'popt': popt,
                    'pcov': pcov
                }
            else:
                fitting_results = {
                    'fit_x': fit_xs,
                    'fit_y': fit_ys,
                    'popt': popt,
                    'pcov': pcov
                }

        return raw_data, fitting_results, num_states

    def get_fisher_info(self, xs,ys):

        fisher_info = np.zeros_like(xs[0], dtype=float)

        for detunings, pops in zip(xs, ys):
            dp_ddetuning = np.gradient(pops, detunings)
            floor = 1e-15
            fisher_info += (dp_ddetuning ** 2) / np.clip(pops, floor, 1)
        return fisher_info



    def fit_with_alpha_multistart(self, detunings, population_vals, fit_funcs, num_states, base_p0, bounds,
                                  fit_x):
        y_data = np.concatenate([
            population_vals[:, state]
            for state in range(num_states)
        ])

        def combined_fit_func(x, d, alpha, t, t0, delta_0, phi):
            return np.concatenate([
                fit_funcs[state](x, d, alpha, t, t0, delta_0, phi)
                for state in range(num_states)
            ])

        alpha_guesses = [0.02, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
        phi_guesses = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]

        best = None

        for alpha0 in alpha_guesses:
            for phi0 in phi_guesses:
                p0 = list(base_p0)
                p0[1] = alpha0  # assuming parameter order: d, alpha, t, t0, delta_0, phi
                p0[5] = phi0

                try:
                    popt, pcov = curve_fit(
                        combined_fit_func,
                        detunings,
                        y_data,
                        p0=p0,
                        bounds=bounds,
                        maxfev=20000,
                    )

                    residual = y_data - combined_fit_func(detunings, *popt)
                    rss = np.sum(residual ** 2)

                    fit_y_list = [
                        fit_funcs[state](fit_x, *popt)
                        for state in range(num_states)
                    ]

                    if best is None or rss < best["rss"]:
                        best = {
                            "popt": popt,
                            "pcov": pcov,
                            "rss": rss,
                            "alpha0": alpha0,
                            'fit_y_list': fit_y_list,
                        }

                except Exception as e:
                    pass

        if best is None:
            fit_y_list = [[None] * len(fit_x) for _ in range(num_states)]
            best = {
                "popt": None,
                "pcov": None,
                "rss": None,
                "alpha0": None,
                'fit_y_list': fit_y_list,
            }

            print("All alpha initial guesses failed")

        return best


    def _make_legend_labels(self, num_states):
        num_ions = num_states - 1
        labels = ["d" * (num_ions - i) + "b" * i for i in range(num_states)]
        labels.append('FI')
        return labels

    def _process_data(self,
                         data_tmp,
                         detuning_idx,
                         secular_idx = 0,
                         use_secular = False):
        # get results
        if use_secular:
            secular_arr = ftw_to_frequency_khz(array(data_tmp[:, secular_idx]))
        else:
            secular_arr = np.zeros(len(data_tmp))

        counts_arr = array(data_tmp[:, 1])
        detuning_arr = ftw_to_frequency_khz(array(data_tmp[:, detuning_idx]))


        secular_vals = np.unique(secular_arr)
        detuning_vals = np.unique(detuning_arr)

        # calculate fluorescence detection threshold
        threshold_list = np.sort(findThresholdScikit(counts_arr))
        num_states = len(threshold_list) + 1

        data_storer = {}

        for sec_idx, sec in enumerate(secular_vals):
            mask = secular_arr == sec

            detuning_subset = detuning_arr[mask]
            counts_subset = counts_arr[mask]

            count_states = np.digitize(counts_subset, threshold_list)

            # group all shots by identical detuning
            detunings, detuning_idx = np.unique(detuning_subset, return_inverse=True)
            population_vals = np.zeros((len(np.unique(detunings)), num_states))
            population_err = np.zeros((len(np.unique(detunings)), num_states))

            for det_idx in range(len(detunings)):
                det_mask = det_idx == detuning_idx
                for state in range(num_states):
                    counts_masked = count_states[det_mask] == state
                    population_vals[det_idx, state] = np.mean(counts_masked)
                    population_err[det_idx, state] = np.std(counts_masked) / np.sqrt(len(counts_masked))


            data_storer[sec_idx] = {
                "secular": sec,
                "detuning": detunings,
                "populations": population_vals,
                "populations_err": population_err
            }

        return data_storer, num_states


    def _get_fit_funcs(self, num_states):
        if num_states == 2:
            fit_funcs = get_single_ion_cat_lineshape()
        elif num_states == 3 and not self.enable_ms_gate:
            fit_funcs = get_unentangled_two_ion_cat_lineshape()
        elif num_states == 3 and self.enable_ms_gate:
            fit_funcs = get_entangled_two_ion_cat_lineshape()
        else:
            raise NotImplementedError
        return fit_funcs