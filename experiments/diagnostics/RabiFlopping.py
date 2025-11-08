from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from sipyco import pyon
from numpy import array, linspace, zeros, int64, where, mean, pi

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, QubitPulseShape, RescueIon, NoOperation, SidebandCoolContinuousRAM
)


class RabiFlopping(LAXExperiment, Experiment):
    """
    Experiment: Rabi Flopping
    Characterize coherent evolution by flopping on a narrow-line transition.
    """
    name = 'Rabi Flopping'
    kernel_invariants = {
        # hardware parmeters
        'freq_rabiflop_ftw', 'ampl_qubit_asf', 'att_readout_mu', 'time_rabiflop_mu_list',

        # subsequences
        'initialize_subsequence', 'doppler_subsequence', 'sbc_subsequence', 'pulseshape_subsequence',
        'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_readout', 'profile_729_SBC',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=80, precision=0, step=1, min=1, max=10000))
        self.setattr_argument("enable_linetrigger", BooleanValue(default=False),
                              tooltip="Trigger the beginning of each shot from the AC line.")
        self.setattr_argument("cooling_type", EnumerationValue(["Doppler", "SBC"], default="SBC"),
                              tooltip="Select the cooling type to use.")

        # rabi flopping arguments
        self.setattr_argument("time_rabi_us_list", Scannable(
                                                    default=[
                                                        RangeScan(0.01, 20, 80, randomize=True),
                                                        ExplicitScan([6.05]),
                                                        CenterScan(3.05, 5., 0.1, randomize=True),
                                                    ],
                                                    global_min=0.01, global_max=100000, global_step=1,
                                                    unit="us", scale=1, precision=3
                                                ), group=self.name)
        self.setattr_argument("freq_rabiflop_mhz", NumberValue(default=101.0968, precision=6, step=1, min=50., max=400., scale=1., unit='MHz'),
                              group=self.name)
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50, precision=3, step=5, min=1, max=50, scale=1., unit='%'),
                              group=self.name)
        self.setattr_argument("att_readout_db", NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, scale=1., unit='dB'),
                              group=self.name)
        self.setattr_argument("equalize_delays", BooleanValue(default=False),
                              group=self.name,
                              tooltip="Adds a dummy delay to each rabi pulse to equalize the time taken for each shot.")
        self.setattr_argument("enable_pulseshaping", BooleanValue(default=False),
                              group=self.name,
                              tooltip="Shape the rabiflop pulse to reduce spectral leakage. "
                                      "Uses a Hann (sine-squared) envelope.")

        # allocate relevant beam profiles
        self.profile_729_readout =  0
        self.profile_729_SBC =      1

        # prepare sequences
        self.sbc_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.pulseshape_subsequence =   QubitPulseShape(
            self, ram_profile=self.profile_729_readout, ram_addr_start=502,
            num_samples=200, ampl_max_pct=self.ampl_qubit_pct, pulse_shape="sine_squared"
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.doppler_subsequence =      NoOperation(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

        # get devices
        self.setattr_device('qubit')
        self.setattr_device('trigger_line')

    def prepare_experiment(self):
        """
        Prepare values for speedy evaluation.
        """
        # choose correct cooling subsequence
        if self.cooling_type == "Doppler":  self.cooling_subsequence = self.doppler_subsequence
        elif self.cooling_type == "SBC":    self.cooling_subsequence = self.sbc_subsequence

        # convert input arguments to machine units
        self.freq_rabiflop_ftw =    self.qubit.frequency_to_ftw(self.freq_rabiflop_mhz * MHz)
        self.ampl_qubit_asf =       self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_readout_mu =       att_to_mu(self.att_readout_db * dB)

        # create timing list such that all shots have same length
        max_time_us = max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list = array([
            [self.core.seconds_to_mu((max_time_us - time_us) * us), self.core.seconds_to_mu(time_us * us)]
            for time_us in self.time_rabi_us_list
        ])
        # turn off delay equalization based on input
        if not self.equalize_delays: self.time_rabiflop_mu_list[:, 0] = int64(8)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # set qubit readout waveform
        if self.enable_pulseshaping:
            self.qubit.set_ftw(self.freq_rabiflop_ftw)
        else:
            self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.ampl_qubit_asf,
                              profile=self.profile_729_readout, phase_mode=PHASE_MODE_CONTINUOUS)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        time_rabi_actual_mu = -1 # holder variable to support variable pulse shaping

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            # sweep rabi flopping times
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:

                '''
                PREPARE SHOT
                '''
                # configure qubit pulse
                self.core.break_realtime()
                if self.enable_pulseshaping:
                    time_rabi_actual_mu = self.pulseshape_subsequence.configure(time_rabi_pair_mu[1])
                    delay_mu(50000)
                else:
                    time_rabi_actual_mu = time_rabi_pair_mu[1]

                # wait for linetrigger
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, self.trigger_line.time_holdoff_mu)


                '''
                RUN SHOT
                '''
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.cooling_subsequence.run_dma()

                # prepare qubit beam for readout
                self.qubit.set_profile(self.profile_729_readout)
                self.qubit.set_att_mu(self.att_readout_mu)

                # add delay to equalize shot time
                delay_mu(time_rabi_pair_mu[0])
                # rabi flop & read out
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.run()
                else:
                    self.qubit.on()
                    delay_mu(time_rabi_actual_mu)
                    self.qubit.off()

                # do readout & clean up loop
                self.readout_subsequence.run_dma()
                self.rescue_subsequence.resuscitate()
                self.initialize_subsequence.slack_rescue()
                counts = self.readout_subsequence.fetch_count()

                # retrieve results & store in dataset
                self.rescue_subsequence.detect_death(counts)
                self.update_results(time_rabi_actual_mu, counts)

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        Fit rabi flopping data with an exponentially damped sine curve.
        """
        # get results
        results_tmp = array(self.results)
        probability_vals = zeros(len(results_tmp))
        counts_arr = array(results_tmp[:, 1])
        # convert x-axis (time) from machine units to seconds
        results_tmp[:, 0] = array([self.core.mu_to_seconds(time_mu) for time_mu in results_tmp[:, 0]])

        # calculate fluorescence detection threshold
        threshold_list = findThresholdScikit(results_tmp[:, 1])
        for threshold_val in threshold_list:
            probability_vals[where(counts_arr > threshold_val)] += 1.
        # normalize probabilities and convert from D-state probability to S-state probability
        results_tmp[:, 1] = 1. - probability_vals / len(threshold_list)
        # process dataset into x, y, with y being averaged probability
        results_tmp = groupBy(results_tmp, column_num=0, reduce_func=mean)
        results_tmp = array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()

        # fit rabi flopping using damped harmonic oscillator
        results_plotting_x, results_plotting_y = array(results_tmp).transpose()
        fit_x = linspace(min(results_plotting_x), max(results_plotting_x), 1000)
        try:
            fitter = fitDampedOscillator()
            fit_params, fit_err = fitter.fit(results_tmp)
            fit_y = fitter.fit_func(fit_x, *fit_params)
            # todo: use fit parameters to attempt to fit roos eqn(A.5)
            # todo: note: we fit using roos' eqn(A.5) instead of eqn(A.3) for simplicity

            # process fit parameters to give values of interest
            fit_period_us = (2 * pi * 1.e6) / fit_params[2]
            fit_period_err_us = fit_period_us * (fit_err[2] / fit_params[2])
            # todo: extract phonon number from fit

            # save results to hdf5 as a dataset
            self.set_dataset('fit_params', fit_params)
            self.set_dataset('fit_err', fit_err)

            # save results to dataset manager for dynamic experiments
            res_dj = [[fit_period_us, fit_period_err_us], [fit_params, fit_err]]
            self.set_dataset('temp.rabiflopping.results', res_dj, broadcast=True, persist=False, archive=False)
            self.set_dataset('temp.rabiflopping.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

            # print out fitted parameters
            print("\tResults - Rabi Flopping:")
            print("\t\tPeriod (us):\t{:.2f} +/- {:.2f}".format(fit_period_us, fit_period_err_us))

        except Exception as e:
            print("\tUnable to Find Fit for Rabi Flopping")
            fit_y = [None]*len(fit_x)

        # format dictionary for applet plotting
        plotting_results = {'x': results_plotting_x * 1e6,
                            'y': results_plotting_y,
                            'fit_x': fit_x * 1e6,
                            'fit_y': fit_y,
                            'subplot_titles': f'Rabi Flopping',
                            'subplot_x_labels': 'Time (us)',
                            'subplot_y_labels': 'S State Population',
                            'rid': self.scheduler.rid,}
        self.set_dataset('temp.plotting.results_rabi_flopping', pyon.encode(plotting_results), broadcast=True)

        # create applet
        self.ccb.issue("create_applet", f"Data Plotting",
                       '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_rabi_flopping'
                       ' --num-subplots 1',
                       group = ['plotting', 'diagnostics'])
        return results_tmp

