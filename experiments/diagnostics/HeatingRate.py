import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.SidebandCooling as SidebandCooling
from sipyco import pyon


class HeatingRate(SidebandCooling.SidebandCooling):
    """
    Experiment: Heating Rate

    Measures the heating rate by doing sideband cooling, then waiting
    a variable amount of time before readout.
    """
    name = 'Heating Rate'
    kernel_invariants = {
        'freq_sideband_readout_ftw_list', 'time_heating_rate_mu_list',
        'config_heating_rate_list'
    }

    def build_experiment(self):
        # heating rate wait times
        self.setattr_argument("time_heating_rate_ms_list", PYONValue([1,2,3]))
        self.setattr_device('ccb')

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list

        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list = np.array([self.core.seconds_to_mu(time_ms * ms)
                                                   for time_ms in self.time_heating_rate_ms_list],
                                                  dtype=np.int64)

        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_heating_rate_list = np.stack(np.meshgrid(self.time_heating_rate_mu_list,
                                                             self.freq_sideband_readout_ftw_list),
                                                 -1).reshape(-1, 2)
        self.config_heating_rate_list = np.array(self.config_heating_rate_list, dtype=np.int64)
        np.random.shuffle(self.config_heating_rate_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_heating_rate_list),
                3)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            for config_vals in self.config_heating_rate_list:
                # extract values from config list
                time_heating_delay_mu = config_vals[0]
                freq_readout_ftw = np.int32(config_vals[1])
                self.core.break_realtime()

                # set frequency for readout
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # wait time to measure heating rate
                delay_mu(time_heating_delay_mu)

                # sideband readout
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # get results & update dataset
                self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
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
        Fit resultant spectra with a sinc profile to extract n,
        then fit a line to extract heating rate
        """
        # create data structures for processing
        results_tmp = np.array(self.results)
        probability_vals = np.zeros(len(results_tmp))
        counts_arr = np.array(results_tmp[:, 1])
        time_readout_us = self.sidebandreadout_subsequence.time_sideband_readout_us
        num_subplots = len(self.time_heating_rate_ms_list) + 1
        # tmp remove
        results_yzde = groupBy(results_tmp, column_num=2)
        self._tmp_data = results_yzde.copy()
        # tmp remove

        # convert column 0 (frequency) from frequency tuning word (FTW) to MHz (in absolute units),
        # and convert column 2 (time) from machine units to seconds
        results_tmp *= np.array([2.e3 / 0xFFFFFFFF, 1., 1.e-9])

        # calculate fluorescence detection threshold
        threshold_list = findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] = 1. - probability_vals / len(threshold_list)
        # process dataset into x, y, with y being averaged probability
        results_tmp = groupBy(results_tmp, column_num=2)

        # batch process sideband cooling data for each heating time
        heating_rate_data = np.array([[heat_time, *(self._extract_phonon(dataset, time_readout_us))]
                                      for heat_time, dataset in results_tmp.items()])


        # fit line to results
        line_fitter = fitLine()
        # fit_params_heating_rate = np.array(line_fitter.fit(heating_rate_data[:, :2], bounds=((0., -np.inf), (1., np.inf))))
        fit_params_heating_rate = np.array(line_fitter.fit(heating_rate_data[:, :2]))
        # save results to hdf5 as a dataset
        self.set_dataset('processed_heating_rate_data', heating_rate_data)
        self.set_dataset('fit_params_heating_rate', fit_params_heating_rate)

        # print out fitted parameters
        print("\tResults - Heating Rate:")
        print("\t---------------------")
        for heat_time_s, phonon_num, phonon_err in heating_rate_data:
            print("\t\t{:.1f}\tms:\t{:.2f} +/- {:.2f}".format(heat_time_s * 1.e3, phonon_num, phonon_err))
        print("\t---------------------")
        print("\t\tSlope:\t\t{:.3f} quanta/s".format(fit_params_heating_rate[1]))
        print("\t\tIntercept:\t{:.3f} quanta".format(fit_params_heating_rate[0]))

        plotting_results_x = []
        plotting_results_y = []
        fit_x = []
        fit_y = []


        dataset_length = 0
        for heat_time, dataset in results_tmp.items():
            data_rsb, data_bsb = np.array(self._extract_populations(dataset, time_readout_us))
            data_rsb_x, data_rsb_y = data_rsb.transpose()
            data_bsb_x, data_bsb_y = data_bsb.transpose()
            plotting_results_x.append(data_rsb_x/2)
            plotting_results_y.append(data_rsb_y)
            plotting_results_x.append(data_bsb_x / 2)
            plotting_results_y.append(data_bsb_y)
        # plotting_results_x.append(list(heating_rate_data[:, 0]))
        # plotting_results_y.append(list(heating_rate_data[:, 1]))
        heating_rate_fit_x = np.linspace(np.min(heating_rate_data[:, 0]), np.max(heating_rate_data[:, 0]),
                                         len(heating_rate_data[:, 0])*10)
        fit_y.append(fit_params_heating_rate[0]*heating_rate_fit_x + fit_params_heating_rate[1])\

        plotting_results = {'x': plotting_results_x,
                            'y': plotting_results_y,
                            'rid': self.scheduler.rid,
                            }

        self.set_dataset('temp.plotting.results_heating_rate', pyon.encode(plotting_results), broadcast=True)

        num_subplots = 6
        ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_heating_rate'
        ccb_command += f' --num-subplots {num_subplots}'
        self.ccb.issue("create_applet", f"Heating Rate",
                       ccb_command,
                       group = 'plotting.diagnostics')

