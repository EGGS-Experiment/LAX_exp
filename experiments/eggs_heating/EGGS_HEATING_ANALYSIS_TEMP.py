from artiq.experiment import *
import artiq
import time
import numpy as np
import matplotlib.pyplot as plt
import h5py
from analysis.processing import *
from analysis.plotting import *
from LAX_exp.analysis import *
from sipyco import pyon
from sipyco.sync_struct import Subscriber
from sipyco.asyncio_tools import atexit_register_coroutine


class AnalyzeEggsHeating(EnvExperiment):
    name = 'EGGS HEATING ANALYSIS TEMP'


    def build(self):
        self.scheduler = self.get_device("scheduler")

    def run(self):
        pass
    def analyze(self):

        # self.analyze_scans()

    def analyze_linewidth(self):

    def analyze_scans(self):

        file_paths = [r'\\eric.physics.ucla.edu\groups\motion\Data\2024-04\2024-04-19\000056795-EGGS Heating.h5']

        # file_paths = [r'\\eric.physics.ucla.edu\groups\motion\Data\2024-04\2024-04-24\000057196-EGGS Heating.h5']
        # file_paths = [r'\\eric.physics.ucla.edu\groups\motion\Data\2024-04\2024-04-19\000056776-EGGS Heating.h5']
        # file_paths = [r'\\eric.physics.ucla.edu\\groups\\motion\\Data\\2024-05\2024-05-03\\000058838-EGGS Heating.h5']

        for file_path in file_paths:
            dataset = []
            with (h5py.File(file_path) as f):
                dataset = np.array(f['datasets']['results'])
                reps = f['arguments'].attrs['repetitions']

                try:
                    sub_reps = f['arguments'].attrs['sub_repetitions']
                except Exception as e:
                    sub_reps = 1

            ## determine which scan we are running
            if len(np.unique(dataset[:, 2])) > 1:  # carrier sweep
                sorting_col_num = 2

            elif len(np.unique(dataset[:, 3])) > 1:  # sideband sweep
                sorting_col_num = 3
            else:  # sideband readout
                sorting_col_num = 0

            ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
                                                                                           1, 0, reps, sub_reps)

            # todo: automatically determine whether to use coherent or squeezed state for conversion
            phonons = convert_ratios_to_coherent_phonons(ratios)

            if sorting_col_num == 3 or sorting_col_num == 2 or sorting_col_num == 4:
                phonon_err = 0
                if sorting_col_num == 3:
                    fit_params_phonon, fit_err_phonon, fit_data = fitSincGeneric(scanning_freq_MHz, phonons)
                    phonon_max = str(np.round(fit_params_phonon[0], 3))
                else:
                    fit_data = None
                    phonon_max = "No Phonon Number Given for EGGs Heating"
                plot_phonons(scanning_freq_MHz, phonons, fit_data)
                plot_sb_probs(scanning_freq_MHz, ave_rsb, ave_bsb, std_rsb, std_bsb)


            else:
                rsb_freqs_MHz, bsb_freqs_MHz, _ = extract_sidebands_freqs(scanning_freq_MHz)
                s_state_pop_rsb = 1. - ave_rsb
                s_state_pop_bsb = 1. - ave_bsb
                fit_params_rsb, fit_err_rsb, fit_rsb = fitSincGeneric(rsb_freqs_MHz, ave_rsb)
                fit_params_bsb, fit_err_bsb, fit_bsb = fitSincGeneric(bsb_freqs_MHz, ave_bsb)
                phonon_max = str(np.round(fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0]), 3))

                print(fit_params_rsb)

                plot_sidebands(rsb_freqs_MHz, bsb_freqs_MHz, s_state_pop_rsb, s_state_pop_bsb, std_rsb, std_bsb,
                               1. - fit_rsb, 1. - fit_bsb)

            print("\tResults - Sideband Cooling:")
            print("\t\tn:\t{}".format(phonon_max))
        plt.show()

