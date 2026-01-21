from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import ParametricExcite

from random import shuffle
from itertools import groupby
from numpy.linalg import norm
from numpy import array, zeros, mean, argsort, linspace, abs, angle, exp, pi

from EGGS_labrad.config.dc_config import dc_config

# todo: implement adaptive option, and improve use of INITIAL_MODE_VECTORS
# todo: need to create a "running store"/queue for host operations that's NOT used for storage - this will simplify host/kernel interactions
# todo: create a dict for "axis config"
# todo: create function called reset_initial_voltages so we can reset upon error condition
# todo: finish documenting tooltips

# todo: should we put all _host variables as kern_inv? they technically are...


class MicromotionCompensation(LAXExperiment, Experiment):
    """
    Experiment: Micromotion Compensation

    Characterize the micromotion along a mode by applying a parametric excitation on the
    trap RF while scanning shim voltages, then attempt to algorithmically compensate for it.
    """
    name = 'Micromotion Compensation'

    kernel_invariants = {
        # voltages & scan axes
        'dc_config_channeldict', "dc_channel_axes_nums", "dc_channel_axes_names",
        "dc_scan_range_volts_list", "time_dc_synchronize_delay_mu",

        # DDS values
        "ampl_cooling_asf", "freq_cooling_ftw", "att_db_bounds",
        "profile_397_parametric", "profile_dds_parametric",

        # magic numbers
        "OPT_CORR_AMPL_FRAC", "GUESS_CORR_AMPL_GAMMA", "CORR_AMPL_ATT_SLOPE", "CONVERGENCE_VOLTAGE_V",
        "INITIAL_MODE_VECTORS", 'TIME_DC_SYNC_MS', 'OPTIMA_LENGTH_FIT', 'MIN_VOLTAGE_SCAN_RANGE',

        # subsequences
        "parametric_subsequence",
    }

    def build_experiment(self):
        # get DC channel configuration dictionary
        self.dc_config_channeldict = dc_config.channeldict

        # core arguments
        self.setattr_argument("iterations", NumberValue(default=3, precision=0, step=1, min=1, max=10),
                              tooltip="The number of optimization iterations to perform.\n"
                                      "Each iteration consists of a voltage scan on each axes, with an optimum update after each voltage scan.\n"
                                      "If successive optima converge within some target value, then this experiment will stop early.")

        # general configuration
        self.setattr_argument("repetitions_per_voltage", NumberValue(default=6, precision=0, step=1, min=1, max=100),
                              group='configuration',
                              tooltip="The number of repetitions on each mode for a single voltage point.\n"
                                      "All repetitions are executed at once for each point, but are alternated "
                                      "between modes (e.g. 2 reps => mode0, mode1, mode0, mode1).")
        self.setattr_argument("num_steps", NumberValue(default=14, precision=0, step=1, min=5, max=100),
                              group='configuration',
                              tooltip="The number of steps for each voltage scan.")
        self.setattr_argument("adaptive", BooleanValue(default=True), group='configuration',
                              tooltip="Predicts the mode vectors for each scan, and automatically adjusts "
                                      "the attenuation.\n"
                                      "This feature is currently always on (i.e. False will do nothing lol).")

        # modulation - mode #1
        self.setattr_argument("freq_mode_0_khz",    NumberValue(default=1575.11, precision=3, step=10, min=1, max=10000, scale=1., unit="kHz"),
                              group='modulation')
        self.setattr_argument("att_mode_0_db",      NumberValue(default=21., precision=1, step=0.5, min=0, max=31.5, scale=1., unit="dB"),
                              group='modulation',
                              tooltip="The starting attenuation to use for mode 0.\n"
                                      "If the adaptive option is selected, the experiment will automatically adjust "
                                      "the attenuation to maximize signal.")
        # modulation - mode #2
        self.setattr_argument("freq_mode_1_khz",    NumberValue(default=1287.55, precision=3, step=10, min=1, max=10000, scale=1., unit="kHz"),
                              group='modulation')
        self.setattr_argument("att_mode_1_db",      NumberValue(default=25., precision=1, step=0.5, min=0, max=31.5, scale=1., unit="dB"),
                              group='modulation',
                              tooltip="The starting attenuation to use for mode 1.\n"
                                      "If the adaptive option is selected, the experiment will automatically adjust "
                                      "the attenuation to maximize signal.")

        # shim voltages
        self.setattr_argument("dc_channel_axis_0",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='V Shim'),
                              group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_0", PYONValue([45, 85]),
                              group='voltages')
        self.setattr_argument("dc_channel_axis_1",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='H Shim'),
                              group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_1", PYONValue([55, 95]), group='voltages')

        # cooling
        self.setattr_argument("ampl_cooling_pct",   NumberValue(default=17, precision=2, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group='cooling',
                              tooltip="The DDS amplitude to use for the cooling beam. Lower amplitudes reduce background "
                                      "and improve SNR, but will take longer and makes the ion prone to death.")
        self.setattr_argument("freq_cooling_mhz",   NumberValue(default=103, precision=6, step=1, min=1, max=500, scale=1., unit="MHz"),
                              group='cooling',
                              tooltip="The DDS frequency to use for the cooling beam. Frequencies closer to resonance "
                                      "improve the SNR, but makes the ion prone to death.")

        # explicitly specify AD9910 profiles
        self.profile_397_parametric = 0
        self.profile_dds_parametric = 6

        # get relevant devices
        self.setattr_device('pmt')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_parametric')
        self.setattr_device('trap_dc')

        # get relevant subsequences
        self.parametric_subsequence =   ParametricExcite(self)


        '''
        MAGIC NUMBERS
        '''
        self.OPT_CORR_AMPL_FRAC =       0.2     # target corr ampl when adjusting atts
        self.GUESS_CORR_AMPL_GAMMA =    0.15    # used to predict target corr ampl
        self.CORR_AMPL_ATT_SLOPE =      -7.362e-3   # guess for dependence of corr ampl on att
        self.CONVERGENCE_VOLTAGE_V =    0.3     # voltage "distance" criteria for convergence
        self.TIME_DC_SYNC_MS =          988     # holdoff period after we set a voltage to allow it to settle
        self.OPTIMA_LENGTH_FIT =        4       # number of prev optima to fit when guessing global opt
        self.MIN_VOLTAGE_SCAN_RANGE =   2       # minimum voltage range for a scan (to ensure we get sufficient excitation)
        self.INITIAL_MODE_VECTORS =     array([-0.5, 0.75]) # initial guess for mode vectors [RF, EGGS]

    def prepare_experiment(self):
        """
        Prepare and predefine variables for the experiment to reduce timing
        overheads during compilation and execution.
        """
        '''
        VOLTAGES
        '''
        # sanitize voltage input
        for voltage_arr in [self.dc_scan_range_volts_axis_0, self.dc_scan_range_volts_axis_1]:
            # check voltage scan ranges have correct format
            if (type(voltage_arr) is not list) or (len(voltage_arr) != 2):
                raise ValueError("Error: voltage scan ranges have incorrect type. Must be list of length 2.")
            # check voltages are all in valid range
            if any([(voltage < 0) or (voltage > 110) for voltage in voltage_arr]):
                raise ValueError("Error: voltage range is out of bounds.")

        # get DC channel numbers & names
        self.dc_channel_axes_nums =     [self.dc_config_channeldict[self.dc_channel_axis_0]['num'],
                                         self.dc_config_channeldict[self.dc_channel_axis_1]['num']]
        self.dc_channel_axes_names =    [self.dc_channel_axis_0, self.dc_channel_axis_1]
        # store voltage bounds conveniently for programmatic access
        self.dc_scan_range_volts_list = [sorted(self.dc_scan_range_volts_axis_0), sorted(self.dc_scan_range_volts_axis_1)]

        # holdoff period after we set a voltage to allow it to settle
        self.time_dc_synchronize_delay_mu = self.core.seconds_to_mu(self.TIME_DC_SYNC_MS * ms)


        '''
        MODULATION
        '''
        # get attenuation bounds from dataset manager
        max_parametric_att_db = self.get_parameter('att_parametric_max_db', group='beams.att_db', override=False)
        if not (0. <= max_parametric_att_db <= 31.5):
            raise ValueError("Max parametric DDS attenuation set incorrectly in dataset manager: {:}.".format(
                max_parametric_att_db))
        self.att_db_bounds = [max_parametric_att_db, 31.5]

        # convert cooling parameters to machine units
        self.ampl_cooling_asf = self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100.)
        self.freq_cooling_ftw = self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)


        '''
        DATA STRUCTURES
        '''
        # holder to store initial voltages in case of booboo
        # note: this isn't kernel_invariant b/c want to grab in initialize_experiment
        #   instead of prepare so that we take values closest to runtime
        self._host_dc_voltages_initial = [0., 0.]

        # create data structure for storing predicted global optima
        # note: only used for convenience purposes (i.e. for later post-processing)
        self.set_dataset('global_optima', zeros((self.iterations * 2 + 1, 2), dtype=float))
        self.setattr_dataset('global_optima')

        # create data structures to store all PROCESSED sweep results on both modes:
        #   i.e. [(v_opt_ax_0_mode0, v_opt_ax_1_mode0, v_fit_err_mode0, slope_ampl_mag_mode0),
        #           (v_opt_ax_0_mode1, v_opt_ax_1_mode1, v_fit_err_mode1, slope_ampl_mag_mode1)]
        # note: used by almost ALL functions
        self.set_dataset('sweep_results', zeros((self.iterations * 2, 2, 4), dtype=float))
        self.setattr_dataset('sweep_results')

        # _host_sweep_counter: index iterator for ALMOST ALL of these data structures
        #   basically used to keep track of time/number of scans we've done
        # note: used by almost ALL functions (i.e. wherever global_optima or sweep_results datasets are used)
        self._host_sweep_counter =  0

        # store most recent predicted optima (note: use center of scan range for starting voltages)
        # todo: can we get rid of this and just use global_optima???
        self._host_voltage_optima_current = array([mean(self.dc_scan_range_volts_axis_0),
                                                   mean(self.dc_scan_range_volts_axis_1)])

        # create data structures to hold the demodulated signal at each voltage position
        # note: used only by _demodulate counts (to store results) and process_optimum (to get results)
        # note: yes, appending to array is inefficient, but reduces complexity of management
        self._host_demod_holder = [[], []]

        # guess initial mode vectors
        # todo: improve & generalize
        if self.freq_mode_0_khz < self.freq_mode_1_khz: self.INITIAL_MODE_VECTORS = array([-0.5, 0.75]) # [RF, EGGS]
        else:   self.INITIAL_MODE_VECTORS = array([0.75, -0.5]) # [EGGS, RF]

    @property
    def results_shape(self):
        # x2 for 2 voltage axes; x2 for 2 modes
        return (self.iterations * 2 * self.num_steps * self.repetitions_per_voltage * 2,
                5)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # prepare host-side values IN AN RPC (b/c values only sync between host and kernel after kernel completion)
        self._initialize_host_values()

        # get DDS CPLD att values so ARTIQ remembers them
        self.core.break_realtime()
        self.dds_parametric.cpld.get_att_mu()

        # prepare cooling beams
        self.core.break_realtime()
        self.pump.rescue()
        self.pump.on()
        self.repump_cooling.on()
        self.repump_qubit.on()
        # store 397nm waveform for modulation
        self.pump.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf,
                         profile=self.profile_397_parametric,
                         phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(10000)

        # set up DDS for modulation
        self.dds_parametric.set_profile(self.profile_dds_parametric)
        delay_mu(8000)

    @rpc
    def _initialize_host_values(self) -> TNone:
        """
        Initialize host values (e.g. pull values from labrad) in an RPC.
        This is necessary since variables are ONLY synchronized between host and kernel at completion.
        """
        # retrieve initial voltages (store in case of booboo later)
        self._host_dc_voltages_initial[0] = self.trap_dc.voltage_get(self.dc_channel_axes_nums[0])
        self._host_dc_voltages_initial[1] = self.trap_dc.voltage_get(self.dc_channel_axes_nums[1])

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        """
        Main sequence of experiment.
        """
        # initial setup!
        self.predict_optimum()  # predict initial optimum
        mode_config = array([[self.freq_mode_0_khz, self.att_mode_0_db],
                             [self.freq_mode_1_khz, self.att_mode_1_db]])   # create a config of [freq_khz, att_db] for each mode

        # run given number of iterations with interleaved optimization along each voltage axis
        for _iter_num in range(self.iterations):
            for scan_axis in range(len(self.dc_channel_axes_nums)):
                # prepare adapted voltage vectors and attenuations for scan
                voltage_scan_arr = self.prepare_voltage_scan(voltage_axis=scan_axis)
                mode_config = self.prepare_att(mode_config, voltage_scan_arr)
                self.scan_voltage_axis(mode_config,
                                       self.dc_channel_axes_nums[scan_axis],
                                       voltage_scan_arr)

                # process optima to check convergence & store results
                self.process_optimum(voltage_axis=scan_axis)
                if self.predict_optimum(): return   # stop early if we converged on optimum
                self.check_termination()


    '''
    HELPER FUNCTIONS
    '''
    @rpc
    def predict_optimum(self) -> TBool:
        """
        Predicts location of global micromotion optimum using linear fits to measured optima of each mode.
        NOTE: interacts with _host_sweep_counter, sweep_results, _host_voltage_optima_current, global_optima
        :return: whether we have converged on the global optimum.
        """
        print("\tBEGIN SWEEP #{:d}: {:}".format(
            self._host_sweep_counter, self.dc_channel_axes_names[self._host_sweep_counter % len(self.dc_channel_axes_nums)]))
        fitter = fitLineLinear()

        # if more than 2 minima, fit line to extract optimum
        if self._host_sweep_counter >= 2:
            # only use relatively recent results to reduce errors
            idx_min = max(self._host_sweep_counter - self.OPTIMA_LENGTH_FIT, 0)
            idx_max = self._host_sweep_counter
            optima_tmp = self.sweep_results[idx_min:idx_max, :, 0:2]

            # fit a line to the optimum for each mode
            fit_mode_0 = fitter.fit(optima_tmp[:, 0])
            fit_mode_1 = fitter.fit(optima_tmp[:, 1])

            # calculate optima as intersection of the fit lines
            opt_v_axis_0 = - (fit_mode_1[0] - fit_mode_0[0]) / (fit_mode_1[1] - fit_mode_0[1])
            opt_v_axis_1 = (fit_mode_1[1] * fit_mode_0[0] - fit_mode_0[1] * fit_mode_1[0]) / (fit_mode_1[1] - fit_mode_0[1])

        # if only have 1 sweep, set optimum as mean of prev sweep and use initial guess to predict other mode's optimum
        elif self._host_sweep_counter == 1:
            # assume voltage axis 0 was scanned
            m0, m1 = self.INITIAL_MODE_VECTORS
            b0, b1 = self.sweep_results[0, :, 0] - self.INITIAL_MODE_VECTORS * mean(self.sweep_results[0, :, 1])
            opt_v_axis_0 = (m1 * b0 - m0 * b1) / (m1 - m0)
            opt_v_axis_1 = - (b1 - b0) / (m1 - m0)

        # simply use current voltages if no previous data exists
        else:
            opt_v_axis_0 = self._host_dc_voltages_initial[0]
            opt_v_axis_1 = self._host_dc_voltages_initial[1]
        print('\t\tPredicted Opt.: ({:.2f} V, {:.2f} V)'.format(opt_v_axis_0, opt_v_axis_1))

        # ensure voltages are within bounds, then update
        if (opt_v_axis_0 < self.dc_scan_range_volts_list[0][0]) or (opt_v_axis_0 > self.dc_scan_range_volts_list[0][1])\
                or (opt_v_axis_1 < self.dc_scan_range_volts_list[1][0]) or (opt_v_axis_1 > self.dc_scan_range_volts_list[1][1]):
            self.reset_initial_voltages()
            raise ValueError("Error: predicted global optimum outside valid scan range.")
        self.trap_dc.voltage_fast(self.dc_channel_axes_nums[0], opt_v_axis_0)
        self.trap_dc.voltage_fast(self.dc_channel_axes_nums[1], opt_v_axis_1)

        # store results
        opt_v_arr = array([opt_v_axis_0, opt_v_axis_1])
        self._host_voltage_optima_current = opt_v_arr
        # todo: why do we use global_optima here instead of e.g. _host_voltage_optima_current?
        self.mutate_dataset('global_optima', self._host_sweep_counter, opt_v_arr)

        # check successive optima for convergence
        if ((self._host_sweep_counter >= 2) and
                (norm(opt_v_arr - self.global_optima[self._host_sweep_counter - 1]) < self.CONVERGENCE_VOLTAGE_V)):
            print('\tSTOP EARLY - CONVERGED ON VOLTAGE OPTIMUM (dist. < {:.2f}).'.format(self.CONVERGENCE_VOLTAGE_V))
            return True
        return False

    @rpc
    def prepare_voltage_scan(self, voltage_axis: TInt32) -> TArray(TFloat, 1):
        """
        Use existing data to generate the voltage list for the sweep.
        NOTE: interacts with sweep_results, _host_voltage_optima_current, _host_sweep_counter
        :param voltage_axis: the voltage axis to be swept.
        :return: an array of voltages to scan.
        """
        # if we have data, use relative difference between previous optima to set the scan range
        if self._host_sweep_counter >= 1:
            # get previous optima
            # todo: do
            opt_prev_mode0 = self.sweep_results[self._host_sweep_counter - 1, 0, 0:2]
            opt_prev_mode1 = self.sweep_results[self._host_sweep_counter - 1, 1, 0:2]

            # set scan range depending on difference between prev optima to try and hit both optima
            prev_axis_num = (voltage_axis - 1) % len(self.dc_channel_axes_nums)
            # todo: why do we use global_optima here instead of e.g. _host_voltage_optima_current?
            voltage_center = self._host_voltage_optima_current[voltage_axis]
            # note: ensure scan range sufficient to get results
            voltage_range = max(self.MIN_VOLTAGE_SCAN_RANGE,
                                1.4 * abs(opt_prev_mode0[prev_axis_num] - opt_prev_mode1[prev_axis_num]))

        # otherwise, scan about existing center, w/~40% of range about existing center
        else:
            voltage_min_v, voltage_max_v = self.dc_scan_range_volts_list[voltage_axis]
            # todo: why do we use global_optima here instead of e.g. _host_voltage_optima_current?
            voltage_center = self.global_optima[self._host_sweep_counter, 0]
            voltage_range = 0.5 * (voltage_max_v - voltage_min_v)

        # set voltage ranges and ensure within bounds
        voltage_min_v = max(voltage_center - 0.5 * voltage_range, self.dc_scan_range_volts_list[voltage_axis][0])
        voltage_max_v = min(voltage_center + 0.5 * voltage_range, self.dc_scan_range_volts_list[voltage_axis][1])
        print("\t\tScan Range (Axis {:d}, {:}): [{:.2f}, {:.2f}]".format(
            voltage_axis, self.dc_channel_axes_names[voltage_axis], voltage_min_v, voltage_max_v))

        # create voltage scan array
        voltage_scan_arr = linspace(voltage_min_v, voltage_max_v, self.num_steps)
        shuffle(voltage_scan_arr)
        return voltage_scan_arr

    @rpc
    def prepare_att(self, mode_configs_all: TArray(TFloat, 2),
                    voltage_scan_v_arr: TArray(TFloat, 1)) -> TArray(TFloat, 2):
        """
        Adjust attenuations for a given scan range to maximize the detected signal.
        This function aims to keep the correlated amplitude at ~0.2 (i.e. 20% correlation).
        NOTE: interacts with sweep_results
        :param mode_configs_all: config of [freq_khz, att_db] for each mode.
        :param voltage_scan_v_arr: the list of voltages to scan.
        :return: the ideal attenuations to use for each mode.
        """
        # extract starting att list from the mode config holder
        att_list_db = [mode_config[1] for mode_config in mode_configs_all]

        # if we have data, use complex fit parameters to guess the optimal attenuation
        if self._host_sweep_counter >= 1:
            # get range of target voltage scan to calculate ideal complex slope magnitude
            voltage_range = (max(voltage_scan_v_arr) - min(voltage_scan_v_arr)) / 2.
            slope_mag_ideal = self.OPT_CORR_AMPL_FRAC / norm((voltage_range, self.GUESS_CORR_AMPL_GAMMA))
            # get most recent fitted complex slope values
            slope_mag_prev = self.sweep_results[self._host_sweep_counter - 1, :, 3]

            # guess ideal attenuations for both modes
            att_db_change = (slope_mag_ideal - slope_mag_prev) / self.CORR_AMPL_ATT_SLOPE
            att_list_db = [val + att_db_change[i] for i, val in enumerate(att_list_db)]

        # otherwise, simply set attenuations to input values
        else:   att_list_db = [self.att_mode_0_db, self.att_mode_1_db]

        # ensure attenuations within bounds before returning
        att_list_db = [min(self.att_db_bounds[1], max(self.att_db_bounds[0], att_db)) for att_db in att_list_db]
        print('\t\tAttenuations:\n'
              '\t\t  Mode 0 ({:.1f} kHz): {:.1f} dB\n\t\t  Mode 1 ({:.1f} kHz): {:.1f} dB'.format(
            mode_configs_all[0, 0], att_list_db[0], mode_configs_all[1, 0], att_list_db[1]))

        # modify original mode_config with updated attenuations
        mode_configs_all[:, 1] = att_list_db
        return mode_configs_all

    @kernel(flags={"fast-math"})
    def scan_voltage_axis(self, mode_configs_all: TArray(TFloat, 2), dc_channel_num: TInt32,
                          voltage_scan_v_arr: TArray(TFloat, 1)) -> TNone:
        """
        Characterize micromotion along a voltage axis via parametric excitation.
        :param mode_configs_all: config of [freq_khz, att_db] for each mode.
        :param dc_channel_num: the voltage channel number to be scanned.
        :param voltage_scan_v_arr: the list of voltages to scan.
        """
        # scan voltage configurations
        for voltage_v in voltage_scan_v_arr:
            # set DC voltage and synchronize hardware clock with timeline
            self.trap_dc.voltage_fast(dc_channel_num, voltage_v)
            self.core.break_realtime() # add slack to guarantee now_mu() in future
            self.core.wait_until_mu(now_mu())
            delay_mu(self.time_dc_synchronize_delay_mu) # add extra delay for voltages to settle

            # run N reps at target voltage (b/c voltage adjustment dominates overheads)
            # and scan over each mode
            for rep_num in range(self.repetitions_per_voltage):
                for mode_idx in range(len(mode_configs_all)):

                    # get mode configs
                    freq_mode_khz = mode_configs_all[mode_idx, 0]
                    att_mode_db =   mode_configs_all[mode_idx, 1]
                    # convert relevant values to machine units
                    # note: yes, would be good to do in advance, but conversion is bothersome and requires keeping
                    #   more global variables, and we're not timing-critical here
                    freq_mode_ftw = self.dds_parametric.frequency_to_ftw(freq_mode_khz * kHz)
                    att_mode_mu =   self.dds_parametric.cpld.att_to_mu(att_mode_db * dB)

                    # prepare hardware for parametric demodulation
                    self.core.break_realtime()
                    self.dds_parametric.set_att_mu(att_mode_mu)
                    self.dds_parametric.set_mu(freq_mode_ftw,
                                               asf=self.dds_parametric.ampl_modulation_asf,
                                               profile=self.profile_dds_parametric,
                                               phase_mode=PHASE_MODE_CONTINUOUS)
                    self.pump.set_profile(self.profile_397_parametric)
                    delay_mu(25000)

                    # run parametric excitation and get timestamps
                    pmt_timestamp_list = self.parametric_subsequence.run()

                    # clean up and process results (e.g. demodulate counts)
                    self.pump.rescue()  # leave beams on rescue while idle
                    self.pmt.clear_inputs()
                    self._demodulate_counts(mode_idx, freq_mode_khz * kHz, voltage_v, pmt_timestamp_list)
                    self.check_termination()

    @rpc
    def _demodulate_counts(self, mode_idx: TInt32, freq_hz: TFloat, voltage_v: TFloat,
                          timestamp_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Demodulate PMT timestamps with respect to the modulation frequency.
        NOTE: interacts with _host_demod_holder
        :param mode_idx: the mode being interrogated.
        :param freq_hz: the modulation frequency (in Hz).
        :param voltage_v: the current electrode voltage (in volts).
        :param timestamp_mu_list: the list of timestamps (in machine units) to demodulate.
        """
        # convert arguments to human units
        timestamps_s = self.core.mu_to_seconds(array(timestamp_mu_list))
        
        # digitally demodulate counts and convert correlated signal to polar coordinates
        correlated_signal = mean(exp((2.j * pi * freq_hz) * timestamps_s))
        correlated_ampl =   abs(correlated_signal)
        correlated_phase =  angle(correlated_signal)
        # extract count rate in seconds
        # todo: actually store lmao
        count_rate_hz =     len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # store results in global parent dataset for post-processing
        res_arr = [freq_hz, voltage_v, correlated_ampl, correlated_phase]
        self.update_results(*res_arr, count_rate_hz)
        self._host_demod_holder[mode_idx].append(res_arr) # store demodulated counts in holder for optimum extraction

    @rpc
    def process_optimum(self, voltage_axis: TInt32) -> TNone:
        """
        Extract voltage optimum for each mode via complex linear fit on demodulated counts.
        This function simplifies host-kernel interactions since raw data is store on host-side.
        This function should only be called at the end of a voltage scan since it updates internal
            counters/iterators (i.e. this function defines the end of a voltage scan).
        NOTE: interacts with _host_demod_holder, _host_sweep_counter, _host_voltage_optima_current
        # todo: check results to ensure max corr ampl is above a given threshold
        :param voltage_axis: the voltage axis that was swept.
        """
        opt_res_store = []  # store processed optima

        # process optima separately for each mode
        for mode_idx, demod_results in enumerate(self._host_demod_holder):
            # sort & group array by voltage (column 1), then average demodulated signal
            # note: sort is necessary here since groupby requires key/values to be contiguous
            res_arr = array(demod_results)
            res_arr_tmp = res_arr[argsort(res_arr[:, 1])]
            res_arr_tmp = array([
                mean(array([val for val in group]), axis=0)
                for key, group in groupby(res_arr_tmp, lambda arr: arr[1])
            ])

            # format results into a 2D array with complex type for complex linear fitting
            results_tmp = array([
                res_arr_tmp[:, 1],
                res_arr_tmp[:, 2] * exp(1.j * res_arr_tmp[:, 3])
            ], dtype='complex128').transpose()
            # sort results by voltage (again)
            results_tmp = results_tmp[argsort(results_tmp[:, 1], axis=0)]

            # extract complex slope of correlated ampl. vs voltage
            b_fit, m_fit = complexLinearFit(results_tmp)
            slope_ampl_mag = norm((m_fit[0], m_fit[1]))

            # extract minimum mode voltage
            opt_voltage_v, opt_voltage_err = complexLinearFitMinimize(results_tmp)
            print("\t\tMode {:d} Opt. (Axis {:d}, {:}): {:.2f} +/- {:.2f} V".format(
                mode_idx, voltage_axis, self.dc_channel_axes_names[voltage_axis],
                opt_voltage_v, opt_voltage_err)
            )

            # check optima for errors
            if ((opt_voltage_v < self.dc_scan_range_volts_list[voltage_axis][0]) or
                    (opt_voltage_v > self.dc_scan_range_volts_list[voltage_axis][1])):
                self.reset_initial_voltages()
                raise ValueError("Error: Mode {:d} voltage out of range: {:f} V.".format(
                    mode_idx, opt_voltage_v))
            elif abs(opt_voltage_err) > 2.0:
                self.reset_initial_voltages()
                raise ValueError("Error: Mode {:d} optimum uncertainty exceeds bounds: {:f} V.".format(
                    mode_idx, opt_voltage_err))

            # get full optima vector and store
            # todo: why do we use global_optima here instead of e.g. _host_voltage_optima_current? or sweep_results???
            opt_v_vector = self._host_voltage_optima_current
            opt_v_vector[voltage_axis] = opt_voltage_v  # only update the voltage axis we just scanned
            opt_res_store.append([*opt_v_vector, opt_voltage_err, slope_ampl_mag])

        # store results, clear holders, and update iterators
        self.mutate_dataset('sweep_results', self._host_sweep_counter, opt_res_store)
        self._host_demod_holder = [[], []]
        self._host_sweep_counter += 1

    @rpc
    def reset_initial_voltages(self) -> TNone:
        """
        Reset shim voltages to initial values.
        Intended to be called in case of error (e.g. voltage optima outside of valid/safe range) to
            return the experiment to a reasonably happy state.
        This function is not async to ensure that it completes immediately.
        """
        self.trap_dc.voltage_set(self.dc_channel_axes_nums[0], self._host_dc_voltages_initial[0])
        self.trap_dc.voltage_set(self.dc_channel_axes_nums[1], self._host_dc_voltages_initial[1])

