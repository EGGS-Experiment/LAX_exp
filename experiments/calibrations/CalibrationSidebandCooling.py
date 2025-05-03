import numpy as np
from sipyco import pyon
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon,
    SidebandCoolContinuousRAM, SidebandCoolPulsed, SidebandReadout
)
# todo: sweep spinpol pct


class CalibrationSidebandCooling(LAXExperiment, Experiment):
    """
    Calibration: Sideband Cooling

    Sideband cooling but with scannable parameters for optimization.
    """
    name = 'Calibration Sideband Cooling'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_readout', 'profile_729_SBC', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("mode_target",    NumberValue(default=1, precision=0, step=1, min=1, max=10))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_readout = 0
        self.profile_729_SBC = 1

        # base SBC
        self.setattr_argument("sideband_cooling_config_list", PYONValue({100.7555: [26., 5.], 100.455: [37., 5.], 100.315: [37., 5.]}),
                              tooltip="{freq_mode_mhz: [sbc_mode_pct_per_cycle, ampl_quench_mode_pct]}", group='SBC.base')
        self.setattr_argument("sideband_cycles_continuous", NumberValue(default=10, precision=0, step=1, min=1, max=10000),
                              tooltip="number of times to loop over the SBC configuration sequence", group='SBC.base')
        self.setattr_argument("time_per_spinpol_us",    NumberValue(default=600, precision=3, step=1, min=0.01, max=100000),
                              tooltip="time between spin polarization pulses (in us)", group='SBC.base')

        # SBC parameter scanning
        self.setattr_argument("freq_sbc_scan_khz_list",  Scannable(
                                                        default=[
                                                            CenterScan(0, 20, 20, randomize=True),
                                                            ExplicitScan([0.]),
                                                            RangeScan(-10, 10, 20, randomize=True),
                                                        ],
                                                        global_min=-1000, global_max=1000, global_step=1.,
                                                        unit="kHz", scale=1, precision=3
                                                    ), group="SBC.sweep")
        self.setattr_argument("time_sbc_us_list",   Scannable(
                                                        default=[
                                                            RangeScan(100, 2000, 100, randomize=True),
                                                            ExplicitScan([6.05]),
                                                            CenterScan(3.05, 5., 0.1, randomize=True),
                                                        ],
                                                        global_min=50, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="SBC.sweep")
        self.setattr_argument("ampl_quench_pct_list",   Scannable(
                                                            default=[
                                                                RangeScan(1., 5., 15., randomize=True),
                                                                ExplicitScan([3.5]),
                                                                CenterScan(3.5, 4., 0.2, randomize=True),
                                                            ],
                                                            global_min=0.01, global_max=50., global_step=1,
                                                            unit="%", scale=1, precision=3
                                                        ), group="SBC.sweep")

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_readout)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # get relevant devices
        self.setattr_device('qubit')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        # todo:
        pass

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # scan over sideband readout frequencies
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # todo

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_readout)
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

