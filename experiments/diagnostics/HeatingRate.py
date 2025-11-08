import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import QubitRAP
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
        # base configs
        'freq_sideband_readout_ftw_list', 'config_experiment_list',

        # RAP configs
        'profile_729_RAP', 'rap_subsequence', 'enable_RAP',
        'att_rap_mu', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',
    }

    def build_experiment(self):
        # heating rate: special arguments
        self.setattr_argument("time_heating_rate_ms_list",  Scannable(
                                                            default=[
                                                                ExplicitScan([1, 10, 20]),
                                                                RangeScan(1, 1000, 200, randomize=True),
                                                            ],
                                                            global_min=0.001, global_max=10000., global_step=10,
                                                            unit="ms", scale=1, precision=3
                                                        ))
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="RAP"))

        # run regular sideband cooling build
        super().build_experiment()
        # extend kernel_invariants w/ SBC parent's kernel_invariants (b/c we redefine for HeatingRate.py)
        kernel_invariants_parent = getattr(super(), "kernel_invariants", set())
        self.kernel_invariants = self.kernel_invariants | kernel_invariants_parent

        # RAP-based readout
        self.profile_729_RAP = 5
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=100.7394, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

        # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=501, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

    def prepare_experiment(self):
        """
        Prepare values for speedy evaluation.
        """
        # run preparations for sideband cooling
        super().prepare_experiment()

        '''CONVERT VALUES TO MACHINE UNITS'''
        time_heating_rate_mu_list = [self.core.seconds_to_mu(time_ms * ms)
                                     for time_ms in self.time_heating_rate_ms_list]

        # prepare RAP values
        self.att_rap_mu = att_to_mu(self.att_rap_db * dB)
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        # configure readout method
        if self.readout_type == 'RAP':
            self.enable_RAP = True
            self.freq_sideband_readout_ftw_list = np.array([self.freq_rap_center_ftw], dtype=np.int32)
        elif self.readout_type == 'Sideband Ratio':
            self.enable_RAP = False
            self.freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        else:
            raise ValueError("Invalid readout type. Must be one of (Sideband Ratio, RAP).")

        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        # create an array of values for the experiment to sweep
        self.config_experiment_list = create_experiment_config(
            time_heating_rate_mu_list, self.freq_sideband_readout_ftw_list,
            shuffle_config=True, config_type=np.int64
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # configure RAP pulse
        # note: initialize RAP here b/c can't call do super().prepare_experiment
        if self.enable_RAP:
            self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                ### PREPARE & CONFIGURE ###
                # extract values from config list
                time_heating_delay_mu = config_vals[0]
                freq_readout_ftw =      np.int32(config_vals[1])

                # set frequency for readout
                self.core.break_realtime()
                if not self.enable_RAP:
                    self.qubit.set_mu(freq_readout_ftw,
                                      asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                      profile=self.profile_729_readout,
                                      phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(10000)


                ### INITIALIZE & DELAY ###
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.sbc_subsequence.run_dma()

                # wait time to measure heating rate
                delay_mu(time_heating_delay_mu)


                ### READOUT ###
                # run readout
                if self.enable_RAP:
                    self.qubit.set_att_mu(self.att_rap_mu)
                    self.rap_subsequence.run_rap(self.time_rap_mu)
                else:
                    self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # clean up loop & update dataset
                self.rescue_subsequence.resuscitate()
                counts = self.readout_subsequence.fetch_count()
                self.initialize_subsequence.slack_rescue()
                self.rescue_subsequence.detect_death(counts)
                self.update_results(freq_readout_ftw, counts, time_heating_delay_mu)

            # rescue ion & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        Fit resultant spectra with a sinc profile to extract n,
        then fit a line to extract heating rate
        """
        # todo: extend analysis to RAP case
        if not self.enable_RAP:
            # create data structures for processing
            results_tmp = np.array(self.results)
            probability_vals = np.zeros(len(results_tmp))
            counts_arr = np.array(results_tmp[:, 1])
            time_readout_us = self.sidebandreadout_subsequence.time_sideband_readout_us
            # tmp remove
            results_yzde = groupBy(results_tmp, column_num=2)
            self._tmp_data = results_yzde.copy()
            # tmp remove
            num_subplots = 2 * len(list(self.time_heating_rate_ms_list)) + 1

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
            try:
                fit_params_heating_rate = np.array(line_fitter.fit(heating_rate_data[:, :2]))
                # print out fitted parameters
                print("\tResults - Heating Rate:")
                print("\t---------------------")
                for heat_time_s, phonon_num, phonon_err in heating_rate_data:
                    print("\t\t{:.1f}\tms:\t{:.2f} +/- {:.2f}".format(heat_time_s * 1.e3, phonon_num, phonon_err))
                print("\t---------------------")
                print("\t\tSlope:\t\t{:.3f} quanta/s".format(fit_params_heating_rate[1]))
                print("\t\tIntercept:\t{:.3f} quanta".format(fit_params_heating_rate[0]))
                self.set_dataset('fit_params_heating_rate', fit_params_heating_rate)
            except Exception as e:
                fit_params_heating_rate = [None, None]
            # save results to hdf5 as a dataset
            self.set_dataset('processed_heating_rate_data', heating_rate_data)

            plotting_results_x = []
            plotting_results_y = []
            fit_x = []
            fit_y = []

            pad_length = int(len(self.freq_sideband_readout_ftw_list)/2 - len(heating_rate_data[:, 0]))
            for heat_time, dataset in results_tmp.items():
                data_rsb, data_bsb = np.array(self._extract_populations(dataset, time_readout_us))
                data_rsb_x, data_rsb_y = data_rsb.transpose()
                data_bsb_x, data_bsb_y = data_bsb.transpose()
                # convert to AOM units
                data_rsb_x /= 2
                data_bsb_x /= 2

                plotting_results_x.append(data_rsb_x)
                plotting_results_y.append(data_rsb_y)
                plotting_results_x.append(data_bsb_x)
                plotting_results_y.append(data_bsb_y)

                ### fit data
                fitter_gauss = fitGaussian()
                fit_x_rsb = np.linspace(np.min(data_rsb_x), np.max(data_rsb_x), len(data_rsb_x)*10)
                fit_x_bsb = np.linspace(np.min(data_bsb_x), np.max(data_bsb_x), len(data_bsb_x)*10)
                try:
                    fit_rsb_params, _ = fitter_gauss.fit(data_rsb)
                    fit_bsb_params, _ = fitter_gauss.fit(data_bsb)
                    fit_y_rsb = fitter_gauss.fit_func(fit_x_rsb, *fit_rsb_params)
                    fit_y_bsb = fitter_gauss.fit_func(fit_x_bsb, *fit_bsb_params)
                    fit_x.append(fit_x_rsb)
                    fit_x.append(fit_x_bsb)
                    fit_y.append(fit_y_rsb)
                    fit_y.append(fit_y_bsb)
                except Exception as e:
                    print("Could not find fit parameters for Sideband")
                    fit_x.append([None]*len(fit_x_rsb))
                    fit_x.append([None]*len(fit_x_bsb))
                    fit_y.append([None]*len(fit_x_rsb))
                    fit_y.append([None]*len(fit_x_bsb))


            plotting_results_x.append(list(heating_rate_data[:, 0]) + [None] * pad_length)
            plotting_results_y.append(list(heating_rate_data[:, 1]) + [None] * pad_length)
            heating_rate_fit_x = np.linspace(np.min(heating_rate_data[:, 0]), np.max(heating_rate_data[:, 0]),
                                             len(heating_rate_data[:, 0]) * 10)
            fit_x.append(heating_rate_fit_x)
            fit_y.append(fit_params_heating_rate[1] * heating_rate_fit_x + fit_params_heating_rate[0])
            ### pad linear fit with Nones
            fit_x[-1] = list(fit_x[-1]) + [None] * (len(fit_x[0])-len(fit_x[-1]))
            fit_y[-1] = list(fit_y[-1]) + [None] * (len(fit_y[0]) - len(fit_y[-1]))

            subplot_titles = np.repeat([f'Wait Time: {round(wait_time,1)}' for wait_time in heating_rate_data[:,0]/ms],2)
            subplot_titles = [str(subplot_title) for subplot_title in subplot_titles]
            plotting_results = {'x': plotting_results_x,
                                'y': plotting_results_y,
                                'fit_x': fit_x,
                                'fit_y': fit_y,
                                'subplot_x_labels': ['AOM Freq (MHz)']*(num_subplots-1)+['Waiting Time (ms)'],
                                'subplot_y_labels': ['D State Population']*(num_subplots-1) + ['Phonon'],
                                'subplot_titles': subplot_titles + ['Summary'],
                                'rid': self.scheduler.rid,
                                'ylims': [[0, 1]]*(num_subplots-1)
                                }

            self.set_dataset('temp.plotting.results_heating_rate', pyon.encode(plotting_results), broadcast=True)
            ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_heating_rate'
            ccb_command += f' --num-subplots {num_subplots}'
            self.ccb.issue("create_applet", f"Data Plotting",
                           ccb_command,
                           group=['plotting', 'diagnostics'])

