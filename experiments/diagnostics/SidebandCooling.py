import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandCoolPulsed, SidebandReadout)


class SidebandCooling(LAXExperiment, Experiment):
    """
    Experiment: Sideband Cooling

    Measures temperature after a given number of RSB pulses.
    """
    name = 'Sideband Cooling'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",        NumberValue(default=10, ndecimals=0, step=1, min=1, max=100000))

        # sideband cooling type
        self.setattr_argument("cooling_type",       EnumerationValue(["Continuous", "Pulsed"], default="Continuous"))

        # get relevant devices
        self.setattr_device('qubit')

        # get subsequences
        self.initialize_subsequence =               InitializeQubit(self)
        self.sidebandcool_pulsed_subsequence =      SidebandCoolPulsed(self)
        self.sidebandcool_continuous_subsequence =  SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =          SidebandReadout(self)
        self.readout_subsequence =                  Readout(self)
        self.rescue_subsequence =                   RescueIon(self)

    def prepare_experiment(self):
        # choose correct cooling subsequence
        if self.cooling_type == "Continuous":       self.sidebandcool_subsequence = self.sidebandcool_continuous_subsequence
        elif self.cooling_type == "Pulsed":         self.sidebandcool_subsequence = self.sidebandcool_pulsed_subsequence

    @property
    def results_shape(self):
        return (self.repetitions * len(self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # scan over sideband readout frequencies
            for freq_ftw in self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # sideband readout
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()


    # ANALYSIS
    def analyze_experiment(self):
        """
        Fit resultant spectrum with a sinc profile.
        """
        # create data structures for processing
        results_tmp =           np.array(self.results)
        probability_vals =      np.zeros(len(results_tmp))
        counts_arr =            np.array(results_tmp[:, 1])
        time_readout_us =       self.sidebandreadout_subsequence.time_sideband_readout_us

        # convert x-axis (frequency) from frequency tuning word (FTW) to MHz (in absolute frequency)
        results_tmp[:, 0] *=    2.e3 / 0xFFFFFFFF


        # calculate fluorescence detection threshold
        threshold_list =        findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] =     1. - probability_vals / len(threshold_list)

        # process dataset into x, y, with y being averaged probability
        results_tmp =           groupBy(results_tmp, column_num=0, reduce_func=np.mean)
        results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()


        # separate spectrum into RSB & BSB and fit using sinc profile
        # guess carrier as mean of highest and lowest frequencies
        guess_carrier_mhz =             (results_tmp[0, 0] + results_tmp[-1, 0]) / 2.
        # split data into RSB and BSB
        split =                         lambda arr, cond: [arr[cond], arr[~cond]]
        results_rsb, results_bsb =      split(results_tmp, results_tmp[:, 0] < guess_carrier_mhz)
        # fit sinc profile
        fit_params_rsb, fit_err_rsb =   fitSinc(results_rsb, time_readout_us)
        fit_params_bsb, fit_err_bsb =   fitSinc(results_bsb, time_readout_us)

        # process fit parameters to give values of interest
        sinc_max =                      lambda a, t: np.sin(np.pi * t * a)**2.
        # phonon_n =                      abs(fit_params_rsb[0]) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))
        phonon_n =                      (abs(sinc_max(fit_params_rsb[0], time_readout_us)) /
                                         (abs(sinc_max(fit_params_bsb[0], time_readout_us)) - abs(sinc_max(fit_params_rsb[0], time_readout_us))))
        phonon_err =                    phonon_n * ((fit_err_rsb[0] / fit_params_rsb[0])**2. +
                                                    (fit_err_rsb[0]**2. + fit_err_bsb[0]**2.) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))**2.
                                                    )**0.5

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
        todo: document
        Arguments:
            ***todo

        Returns:
            ***todo
        """
        # process dataset into x, y, with y being averaged probability
        results_tmp =           groupBy(dataset, column_num=0, reduce_func=np.mean)
        results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()
        time_readout_us =       self.sidebandreadout_subsequence.time_sideband_readout_us



        # separate spectrum into RSB & BSB and fit using sinc profile
        # guess carrier as mean of highest and lowest frequencies
        guess_carrier_mhz =     (results_tmp[0, 0] + results_tmp[-1, 0]) / 2.
        # split data into RSB and BSB
        split =                         lambda arr, cond: [arr[cond], arr[~cond]]
        results_rsb, results_bsb =      split(results_tmp, results_tmp[:, 0] < guess_carrier_mhz)

        # fit sinc profile
        fit_params_rsb, fit_err_rsb =   fitSinc(results_rsb, time_fit_us)
        fit_params_bsb, fit_err_bsb =   fitSinc(results_bsb, time_fit_us)

        # process fit parameters to give values of interest
        sinc_max =                      lambda a, t: np.sin(np.pi*t*a)**2.
        # phonon_n =                      abs(fit_params_rsb[0]) / (abs(fit_params_bsb[0]) - abs(fit_params_rsb[0]))
        phonon_n =                      (abs(sinc_max(fit_params_rsb[0], time_readout_us)) /
                                         (abs(sinc_max(fit_params_bsb[0], time_readout_us)) - abs(sinc_max(fit_params_rsb[0], time_readout_us))))
        phonon_err =                    phonon_n * ((fit_err_rsb[0] / fit_params_rsb[0])**2. +
                                                    (fit_err_rsb[0]**2. + fit_err_bsb[0]**2.) / (fit_params_bsb[0] - fit_params_rsb[0])**2.
                                                    )**0.5
        return np.array([abs(phonon_n), abs(phonon_err)])

