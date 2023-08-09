import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling

# tmp testing
from LAX_exp.analysis import *


class HeatingRate(SidebandCooling.SidebandCooling):
    """
    Experiment: Heating Rate

    Measures the heating rate by doing sideband cooling, then waiting
    a variable amount of time before readout.
    """
    name = 'Heating Rate'


    def build_experiment(self):
        # heating rate wait times
        self.setattr_argument("time_heating_rate_ms_list",                      PYONValue([1, 2, 5, 10, 50]))

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list =                                        np.array([self.core.seconds_to_mu(time_ms * ms)
                                                                                          for time_ms in self.time_heating_rate_ms_list],
                                                                                         dtype=np.int64)

        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_heating_rate_list =                                         np.stack(np.meshgrid(self.time_heating_rate_mu_list, self.freq_readout_ftw_list), -1).reshape(-1, 2)
        self.config_heating_rate_list =                                         np.array(self.config_heating_rate_list, dtype=np.int64)
        np.random.shuffle(self.config_heating_rate_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_heating_rate_list),
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep experiment config: heating time and readout frequency
            for config_vals in self.config_heating_rate_list:

                # extract values from config list
                time_heating_delay_mu =     config_vals[0]
                freq_readout_ftw =          np.int32(config_vals[1])
                self.core.break_realtime()

                # set frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # wait time to measure heating rate
                delay_mu(time_heating_delay_mu)

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)


    # ANALYSIS
    def analyze(self):
        """
        Fit resultant spectra with a sinc profile to extract n,
        then fit a line to extract heating rate
        """
        # create data structures for processing
        results_tmp =           np.array(self.results)
        probability_vals =      np.zeros(len(results_tmp))
        counts_arr =            np.array(results_tmp[:, 1])

        # convert column 0 (frequency) from frequency tuning word (FTW) to MHz,
        # and convert column 2 (time) from machine units to seconds
        results_tmp *=          np.array([1.e3 / 0xFFFFFFFF, 1., 1.e-9])


        # calculate fluorescence detection threshold
        threshold_list =        findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] =     1. - probability_vals / len(threshold_list)
        # process dataset into x, y, with y being averaged probability
        results_tmp =           groupBy(results_tmp, column_num=2)


        # batch process sideband cooling data for each heating time
        heating_rate_data =     np.array([[heat_time, *(self._extract_phonon(dataset, self.time_readout_pipulse_us))]
                                          for heat_time, dataset in results_tmp.items()])
        # fit line to results
        fit_params_heating_rate =       np.array(fitLine(heating_rate_data[:, :2]))

        # save results to hdf5 as a dataset
        self.set_dataset('processed_heating_rate_data', heating_rate_data)
        self.set_dataset('fit_params_heating_rate',     fit_params_heating_rate)

        # print out fitted parameters
        print("\tResults - Heating Rate:")
        for heat_time_ms, phonon_num, phonon_err in heating_rate_data:
            print("\t\t{:.0f} ms:\t{:.2f} +/- {:.2f}".format(heat_time_ms, phonon_num, phonon_err))
        print("\t-------------")
        print("\t\tSlope:\t{:.3f} quanta/s".format(fit_params_heating_rate[0]))
        print("\t\tIntercept:\t{:.3f} quanta".format(fit_params_heating_rate[1]))

    def _extract_phonon(self, dataset, time_fit_us):
        """
        idk
        todo: document
        Arguments:
            ***todo

        Returns:
            ***todo
        """
        # process dataset into x, y, with y being averaged probability
        results_tmp =           groupBy(dataset, column_num=0, reduce_func=np.mean)
        results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()

        # separate spectrum into RSB & BSB and fit using sinc profile
        # guess carrier as mean of highest and lowest frequencies
        guess_carrier_mhz =     (results_tmp[0, 0] + results_tmp[-1, 0]) / 2.
        # split data into RSB and BSB
        def split(arr, cond):
            return [arr[cond], arr[~cond]]
        results_rsb, results_bsb =      split(results_tmp, results_tmp[:, 0] < guess_carrier_mhz)

        # fit sinc profile
        fit_params_rsb, fit_err_rsb =   fitSinc(results_rsb, time_fit_us)
        fit_params_bsb, fit_err_bsb =   fitSinc(results_bsb, time_fit_us)

        # process fit parameters to give values of interest
        phonon_n =                      fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
        phonon_err =                    phonon_n * ((fit_err_rsb[0] / fit_params_rsb[0])**2. +
                                                    (fit_err_rsb[0]**2. + fit_err_bsb[0]**2.) / (fit_params_bsb[0] - fit_params_rsb[0])**2.
                                                    )**0.5
        return np.array([phonon_n, phonon_err])
