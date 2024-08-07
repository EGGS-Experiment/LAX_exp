import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon
from LAX_exp.system.subsequences import NoOperation, SidebandCoolContinuous, SidebandCoolPulsed


class RabiFlopping(LAXExperiment, Experiment):
    """
    Experiment: Rabi Flopping

    Measures ion fluorescence vs 729nm pulse time and frequency.
    """
    name = 'Rabi Flopping'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",            NumberValue(default=30, ndecimals=0, step=1, min=1, max=10000))

        # rabi flopping arguments
        self.setattr_argument("cooling_type",           EnumerationValue(["Doppler", "SBC - Continuous", "SBC - Pulsed"], default="Doppler"))
        self.setattr_argument("time_rabi_us_list",      Scannable(
                                                            default=[
                                                                ExplicitScan([6.05]),
                                                                RangeScan(1, 50, 200, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, ndecimals=5
                                                        ), group=self.name)
        self.setattr_argument("freq_rabiflop_mhz",      NumberValue(default=102.1020, ndecimals=5, step=1, min=1, max=10000), group=self.name)
        self.setattr_argument("att_readout_db",         NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)
        self.setattr_argument("equalize_delays",        BooleanValue(default=True), group=self.name)


        # get devices
        self.setattr_device('qubit')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.doppler_subsequence =                                  NoOperation(self)
        self.sidebandcool_pulsed_subsequence =                      SidebandCoolPulsed(self)
        self.sidebandcool_continuous_subsequence =                  SidebandCoolContinuous(self)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # choose correct cooling subsequence
        if self.cooling_type == "Doppler":                  self.cooling_subsequence =  self.doppler_subsequence
        elif self.cooling_type == "SBC - Continuous":       self.cooling_subsequence =  self.sidebandcool_continuous_subsequence
        elif self.cooling_type == "SBC - Pulsed":           self.cooling_subsequence =  self.sidebandcool_pulsed_subsequence

        # convert input arguments to machine units
        self.freq_rabiflop_ftw =    hz_to_ftw(self.freq_rabiflop_mhz * MHz)
        self.att_readout_mu =       att_to_mu(self.att_readout_db * dB)

        # convert time to machine units
        max_time_us =                   np.max(list(self.time_rabi_us_list))
        # create timing list such that all shots have same length
        self.time_rabiflop_mu_list =    np.array([
                                            [seconds_to_mu((max_time_us - time_us) * us), seconds_to_mu(time_us * us)]
                                            for time_us in self.time_rabi_us_list
                                        ])
        # turn off delay equalization based on input
        if self.equalize_delays:        self.time_rabiflop_mu_list[:, 0] = np.int64(8)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # set qubit readout waveform
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # sweep time
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # run sideband cooling
                self.cooling_subsequence.run_dma()

                # prepare qubit beam for readout
                self.qubit.set_profile(0)
                self.qubit.set_att_mu(self.att_readout_mu)

                # add delay to ensure each shot takes same time
                delay_mu(time_rabi_pair_mu[0])

                # rabi flop
                self.qubit.on()
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.off()

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(time_rabi_pair_mu[1], self.readout_subsequence.fetch_count())
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
        Fit rabi flopping data with an exponentially damped sine curve
        """
        # create data structures for processing
        results_tmp =           np.array(self.results)
        probability_vals =      np.zeros(len(results_tmp))
        counts_arr =            np.array(results_tmp[:, 1])

        # convert x-axis (time) from machine units to seconds
        results_tmp[:, 0] =     np.array([self.core.mu_to_seconds(time_mu) for time_mu in results_tmp[:, 0]])


        # calculate fluorescence detection threshold
        threshold_list =        findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] =     1. - probability_vals / len(threshold_list)

        # process dataset into x, y, with y being averaged probability
        results_tmp =           groupBy(results_tmp, column_num=0, reduce_func=np.mean)
        results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()


        # fit rabi flopping using damped harmonic oscillator
        fit_params, fit_err =   fitDampedOscillator(results_tmp)
        # todo: use fit parameters to attempt to fit roos eqn(A.5)
        # todo: note: we fit using roos' eqn(A.5) instead of eqn(A.3) for simplicity

        # process fit parameters to give values of interest
        fit_period_us =         (2 * np.pi * 1.e6) / fit_params[2]
        fit_period_err_us =     fit_period_us * (fit_err[2] / fit_params[2])
        # todo: extract phonon number from fit

        # save results to hdf5 as a dataset
        self.set_dataset('fit_params',  fit_params)
        self.set_dataset('fit_err',     fit_err)

        # save results to dataset manager for dynamic experiments
        res_dj = [[fit_period_us, fit_period_err_us], [fit_params, fit_err]]
        self.set_dataset('temp.rabiflopping.results', res_dj, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.rabiflopping.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

        # print out fitted parameters
        print("\tResults - Rabi Flopping:")
        print("\t\tPeriod (us):\t{:.2f} +/- {:.2f}".format(fit_period_us, fit_period_err_us))
        return results_tmp
