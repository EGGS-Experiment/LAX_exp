import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon,
    SidebandCoolContinuousRAM, SidebandCoolPulsed, SidebandReadout
)
from sipyco import pyon


class SidebandCooling(LAXExperiment, Experiment):
    """
    Experiment: Sideband Cooling

    Measure ion temperature via sideband ratio following sideband cooling.
    """
    name = 'Sideband Cooling'
    kernel_invariants = {
        'initialize_subsequence', 'sidebandcool_continuous_subsequence', 'sidebandcool_pulsed_subsequence',
        'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))

        # sideband cooling type
        self.setattr_argument("cooling_type",   EnumerationValue(["Continuous", "Pulsed"],
                                                                 default="Continuous"))

        # get relevant devices
        self.setattr_device('qubit')

        # get subsequences
        self.initialize_subsequence =               InitializeQubit(self)
        self.sidebandcool_pulsed_subsequence =      SidebandCoolPulsed(self)
        self.sidebandcool_continuous_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=1, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=500
        )
        self.sidebandreadout_subsequence =          SidebandReadout(self, profile_dds=0)
        self.readout_subsequence =                  Readout(self)
        self.rescue_subsequence =                   RescueIon(self)

    def prepare_experiment(self):
        # choose correct cooling subsequence
        if self.cooling_type == "Continuous":
            self.sidebandcool_subsequence = self.sidebandcool_continuous_subsequence
        elif self.cooling_type == "Pulsed":
            self.sidebandcool_subsequence = self.sidebandcool_pulsed_subsequence

        # shuffle sideband readout frequencies
        np.random.shuffle(self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()
        for trial_num in range(self.repetitions):

            # scan over sideband readout frequencies
            for freq_ftw in self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state & SBC to the ground motional state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # sideband readout
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    # ANALYSIS
    def analyze_experiment(self):
        """
        Fit resultant spectrum with a sinc profile.
        """
        # create data structures for processing
        results_tmp =       np.array(self.results)
        probability_vals =  np.zeros(len(results_tmp))
        counts_arr =        np.array(results_tmp[:, 1])
        time_readout_us =   self.sidebandreadout_subsequence.time_sideband_readout_us
        # convert x-axis (frequency) from frequency tuning word (FTW) to MHz (in absolute frequency)
        results_tmp[:, 0] *= 2.e3 / 0xFFFFFFFF

        # calculate fluorescence detection threshold
        threshold_list = findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] = 1. - probability_vals / len(threshold_list)
        # process dataset into x, y, with y being averaged probability
        results_tmp =   groupBy(results_tmp, column_num=0, reduce_func=np.mean)
        results_tmp =   np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()

        # separate spectrum into RSB & BSB and fit using sinc profile
        # guess carrier as mean of highest and lowest frequencies
        guess_carrier_mhz =             (results_tmp[0, 0] + results_tmp[-1, 0]) / 2.
        # split data into RSB and BSB
        split = lambda arr, cond: [arr[cond], arr[~cond]]
        results_rsb, results_bsb =      split(results_tmp, results_tmp[:, 0] < guess_carrier_mhz)

        # split results into x and y and red and blue sideband
        results_plotting_x_rsb, results_plotting_y_rsb = np.array(results_rsb).transpose()
        results_plotting_x_bsb, results_plotting_y_bsb = np.array(results_bsb).transpose()

        # format arrays for fitting
        fit_x_rsb = np.linspace(np.min(results_plotting_x_rsb), np.max(results_plotting_x_rsb), 1000)
        fit_x_bsb = np.linspace(np.min(results_plotting_x_bsb), np.max(results_plotting_x_bsb), 1000)

        try:
            # attempt to fit sinc profile
            fitter = fitSinc()
            fit_params_rsb, fit_err_rsb = fitter.fit(results_rsb, time_readout_us)
            fit_params_bsb, fit_err_bsb = fitter.fit(results_bsb, time_readout_us)
            fit_y_rsb = fitter.fit_func(fit_x_rsb, *fit_params_rsb)
            fit_y_bsb = fitter.fit_func(fit_x_bsb, *fit_params_bsb)
        except Exception as e:
            # fill array with Nones to avoid fit if fit can't be found
            fit_y_rsb = [None]*len(fit_x_rsb)
            fit_y_bsb = [None]*len(fit_x_bsb)

        # process fit parameters to give values of interest
        sinc_max =  lambda a, t: np.sin(np.pi * t * a)**2.
        # phonon_n =                      abs(fit_params_rsb[0]) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))
        phonon_n =      (abs(sinc_max(fit_params_rsb[0], time_readout_us)) /
                         (abs(sinc_max(fit_params_bsb[0], time_readout_us)) - abs(sinc_max(fit_params_rsb[0], time_readout_us))))
        phonon_err =    phonon_n * ((fit_err_rsb[0] / fit_params_rsb[0])**2. +
                                    (fit_err_rsb[0]**2. + fit_err_bsb[0]**2.) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))**2.
                                    )**0.5

        # format dictionary for applet plotting
        plotting_results = {'x': [results_plotting_x_rsb / 2., results_plotting_x_bsb / 2.],
                            'y': [results_plotting_y_rsb, results_plotting_y_bsb],
                            'fit_x': [fit_x_rsb / 2., fit_x_bsb / 2.],
                            'fit_y': [fit_y_rsb, fit_y_bsb],
                            'subplot_x_labels': 'AOM Freq (MHz)',
                            'subplot_y_labels': 'D State Population',
                            'subplot_titles': ['RSB', 'BSB'],
                            'rid': self.scheduler.rid,
                            'ylims': [[0,1], [0,1]]
                            }

        self.set_dataset('temp.plotting.results_sideband_cooling', pyon.encode(plotting_results), broadcast=True)

        # create applet
        self.ccb.issue("create_applet", f"Sideband Cooling",
                       '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_sideband_cooling'
                       ' --num-subplots 2',
                       group=['plotting', 'diagnostics'])

        # save results to dataset manager for dynamic experiments
        res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
        self.set_dataset('temp.sidebandcooling.results', res_dj, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.sidebandcooling.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
        # save results to hdf5 as a dataset
        self.set_dataset('fit_params_rsb',  fit_params_rsb)
        self.set_dataset('fit_params_bsb',  fit_params_bsb)
        self.set_dataset('fit_err_rsb',     fit_err_rsb)
        self.set_dataset('fit_err_bsb',     fit_err_bsb)
        # print out fitted parameters
        print("\tResults - Sideband Cooling:")
        print("\t\tn: {:.3f} +/- {:.3f}".format(phonon_n, phonon_err))
        print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2., float(fit_err_rsb[1]) / 2.))
        print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2., float(fit_err_bsb[1]) / 2.))
        return results_tmp

    def _extract_phonon(self, dataset, time_fit_us):
        """
        Process a 2D dataset with both rsb and bsb data to extract the phonon number.
        Arguments:
            dataset: array containing results from the experiment
            time_fit_us: time in microseconds experiment was run for

        Returns:
            numpy array containing phonon number and phonon std
        """
        time_readout_us = self.sidebandreadout_subsequence.time_sideband_readout_us
        results_rsb, results_bsb = self._extract_populations(dataset, time_fit_us)
        # fit sinc profile
        fitter = fitSinc()
        fit_params_rsb, fit_err_rsb = fitter.fit(results_rsb, time_fit_us)
        fit_params_bsb, fit_err_bsb = fitter.fit(results_bsb, time_fit_us)

        # process fit parameters to give values of interest
        sinc_max = lambda a, t: np.sin(np.pi*t*a)**2.
        # phonon_n =      abs(fit_params_rsb[0]) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))
        phonon_n =      (abs(sinc_max(fit_params_rsb[0], time_readout_us)) /
                        (abs(sinc_max(fit_params_bsb[0], time_readout_us)) - abs(sinc_max(fit_params_rsb[0], time_readout_us))))
        phonon_err =    phonon_n * ((fit_err_rsb[0] / fit_params_rsb[0])**2. +
                                    (fit_err_rsb[0]**2. + fit_err_bsb[0]**2.) / (fit_params_bsb[0] - fit_params_rsb[0])**2.
                                    )**0.5
        return np.array([abs(phonon_n), abs(phonon_err)])

    def _extract_populations(self, dataset, time_fit_us):
        """
        Extract the population

        Args:
            dataset: array containing results from the experiment
            time_fit_us: time in microseconds experimental trial was run for
        Returns:
            results from red and blue sideband
        """
        # process dataset into x, y, with y being averaged probability
        results_tmp = groupBy(dataset, column_num=0, reduce_func=np.mean)
        results_tmp = np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()

        # separate spectrum into RSB & BSB and fit using sinc profile
        # guess carrier as mean of highest and lowest frequencies
        guess_carrier_mhz = (results_tmp[0, 0] + results_tmp[-1, 0]) / 2.
        # split data into RSB and BSB
        split = lambda arr, cond: [arr[cond], arr[~cond]]
        results_rsb, results_bsb = split(results_tmp, results_tmp[:, 0] < guess_carrier_mhz)

        return results_rsb, results_bsb

