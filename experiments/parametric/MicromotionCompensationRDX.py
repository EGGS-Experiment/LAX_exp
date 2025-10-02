from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import ParametricExcite
import LAX_exp.experiments.parametric.ParametricSweep as ParametricSweep

from random import shuffle
from itertools import groupby
from numpy.linalg import norm
from numpy import array, int32, zeros, mean, argsort, linspace, abs, angle, exp, pi

from EGGS_labrad.config.dc_config import dc_config

# todo: actually make adaptive option usable
# todo: what if we gradient descent within each voltage scan (vs just taking points across the range)?
# todo: improve use of INITIAL_MODE_VECTORS/integrate with adaptive option


class MicromotionCompensation(ParametricSweep.ParametricSweep, Experiment):
    """
    Experiment: Micromotion Compensation

    Characterize the micromotion along a mode by applying a parametric excitation on the
    trap RF while scanning shim voltages, then attempt to algorithmically compensate for it.
    """
    name = 'Micromotion Compensation'

    kernel_invariants = {
        # modulation & modes
        "freq_mode_0_ftw", "freq_mode_1_ftw", "freq_mode_ftw_list", "mode_list_idx",
        "att_mode_0_mu", "att_mode_1_mu", "att_db_bounds_list",

        # voltages
        "dc_channel_axis_0_num", "dc_channel_axis_1_num", "dc_scan_range_volts_list", "time_dc_synchronize_delay_mu",
        "dc_channel_axes_names",

        # beam values
        "ampl_cooling_asf", "freq_cooling_ftw", "profile_397_parametric", "profile_dds_parametric",

        # magic numbers
        "OPT_CORR_AMPL_FRAC", "GUESS_CORR_AMPL_GAMMA", "CORR_AMPL_ATT_SLOPE", "CONVERGENCE_VOLTAGE_V",
        "INITIAL_MODE_VECTORS",

        # subsequences
        "parametric_subsequence",
    }

    def build_experiment(self):
        # get DC channel configuration dictionary
        self.dc_config_channeldict =    dc_config.channeldict
        # explicitly specify AD9910 profiles
        self.profile_397_parametric = 0
        self.profile_dds_parametric = 6

        # core arguments
        self.setattr_argument("iterations", NumberValue(default=3, precision=0, step=1, min=1, max=10))

        # general configuration
        self.setattr_argument("repetitions_per_voltage", NumberValue(default=3, precision=0, step=1, min=1, max=100),
                              group='configuration',
                              tooltip="The number of repetitions on each mode for a single voltage point.\n"
                                      "All repetitions are executed at once for each point, but are alternated "
                                      "between modes (e.g. 2 reps => mode0, mode1, mode0, mode1).")
        self.setattr_argument("num_steps", NumberValue(default=11, precision=0, step=1, min=5, max=100),
                              group='configuration',
                              tooltip="The number of steps for each voltage scan.")
        self.setattr_argument("adaptive", BooleanValue(default=True), group='configuration',
                              tooltip="Predicts the mode vectors for each scan, and automatically adjusts "
                                      "the attenuation. This feature is currently always on (i.e. False will "
                                      "do nothing lol).")

        # modulation - mode #1
        self.setattr_argument("freq_mode_0_khz",    NumberValue(default=1585.22, precision=3, step=10, min=1, max=10000, scale=1., unit="kHz"),
                              group='modulation')
        self.setattr_argument("att_mode_0_db",      NumberValue(default=24., precision=1, step=0.5, min=0, max=31.5, scale=1., unit="dB"),
                              group='modulation',
                              tooltip="The starting attenuation to use for mode 0.\n"
                                      "If the adaptive option is selected, the experiment will automatically adjust "
                                      "the attenuation to maximize signal.")
        # modulation - mode #2
        self.setattr_argument("freq_mode_1_khz",    NumberValue(default=1299.81, precision=3, step=10, min=1, max=10000, scale=1., unit="kHz"),
                              group='modulation')
        self.setattr_argument("att_mode_1_db",      NumberValue(default=27., precision=1, step=0.5, min=0, max=31.5, scale=1., unit="dB"),
                              group='modulation',
                              tooltip="The starting attenuation to use for mode 1.\n"
                                      "If the adaptive option is selected, the experiment will automatically adjust "
                                      "the attenuation to maximize signal.")

        # shim voltages
        self.setattr_argument("dc_channel_axis_0",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='V Shim'),
                              group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_0", PYONValue([50, 90]),
                              group='voltages')
        self.setattr_argument("dc_channel_axis_1",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='H Shim'),
                              group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_1", PYONValue([30, 70]), group='voltages')

        # cooling
        self.setattr_argument("ampl_cooling_pct",   NumberValue(default=18, precision=2, step=5, min=0.01, max=50, scale=1., unit="%"), group='cooling',
                              tooltip="The DDS amplitude to use for the cooling beam. Lower amplitudes reduce background "
                                      "and improve SNR, but will take longer and makes the ion prone to death.")
        self.setattr_argument("freq_cooling_mhz",   NumberValue(default=112, precision=6, step=1, min=1, max=500, scale=1., unit="MHz"), group='cooling',
                              tooltip="The DDS frequency to use for the cooling beam. Frequencies closer to resonance "
                                      "improve the SNR, but makes the ion prone to death.")

        # get relevant devices
        self.setattr_device('pmt')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_parametric')
        self.setattr_device('trap_dc')

        # get relevant subsequences
        self.parametric_subsequence =   ParametricExcite(self)

        '''MAGIC NUMBERS'''
        self.OPT_CORR_AMPL_FRAC =       0.2     # target corr ampl when adjusting atts
        self.GUESS_CORR_AMPL_GAMMA =    0.15    # used to predict target corr ampl
        self.CORR_AMPL_ATT_SLOPE =      -7.362e-3   # guess for dependence of corr ampl on att
        self.CONVERGENCE_VOLTAGE_V =    0.3     # voltage "distance" criteria for convergence
        # self.INITIAL_MODE_VECTORS =     array([-1., 1.]) # initial guess for mode vectors [RF, EGGS]

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
                raise ValueError("InputError: voltage scan ranges have incorrect type.")

            # check voltages are all in valid range
            if any([(voltage < 0) or (voltage > 110) for voltage in voltage_arr]):
                raise ValueError("InputError: voltage range is out of bounds.")

        # get DC channel numbers & names
        self.dc_channel_axis_0_num =    self.dc_config_channeldict[self.dc_channel_axis_0]['num']
        self.dc_channel_axis_1_num =    self.dc_config_channeldict[self.dc_channel_axis_1]['num']
        self.dc_channel_axes_names =    [self.dc_channel_axis_0, self.dc_channel_axis_1]

        # store voltage bounds conveniently for programmatic access
        self.dc_scan_range_volts_list = array([sorted(self.dc_scan_range_volts_axis_0),
                                                  sorted(self.dc_scan_range_volts_axis_1)])

        # holdoff period after we set a voltage to allow it to settle
        self.time_dc_synchronize_delay_mu = self.core.seconds_to_mu(988 * ms)

        '''
        MODULATION
        '''
        # convert modulation (mode #0) parameters to machine units
        self.freq_mode_0_ftw =  self.dds_parametric.frequency_to_ftw(self.freq_mode_0_khz * kHz)
        self.att_mode_0_mu =    att_to_mu(self.att_mode_0_db * dB)

        # convert modulation (mode #1) parameters to machine units
        self.freq_mode_1_ftw =  self.dds_parametric.frequency_to_ftw(self.freq_mode_1_khz * kHz)
        self.att_mode_1_mu =    att_to_mu(self.att_mode_1_db * dB)

        # store modulation values conveniently for programmatic access
        self.freq_mode_ftw_list =   array([self.freq_mode_0_ftw, self.freq_mode_1_ftw], dtype=int32)
        self.mode_list_idx =        list(range(len(self.freq_mode_ftw_list)))
        # get attenuation bounds from dataset
        max_parametric_att_db = self.get_parameter('att_parametric_max_db', group='beams.att_db', override=False)
        if not (0. <= max_parametric_att_db <= 31.5):
            raise ValueError("Max parametric DDS attenuation set incorrectly in dataset manager: {:}.".format(max_parametric_att_db))
        self.att_db_bounds_list = array([max_parametric_att_db, 31.5])

        '''
        COOLING
        '''
        # convert cooling parameters to machine units
        self.ampl_cooling_asf =         self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100.)
        self.freq_cooling_ftw =         self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)

        '''
        DATA STRUCTURES
        '''
        # _host_sweep_counter: index iterator for ALMOST ALL of these data structures
        self._host_sweep_counter =  0

        # create data structure for storing predicted global optima
        # note: only used for convenience purposes (i.e. for later post-processing)
        self.set_dataset('global_optima', zeros((self.iterations * 2 + 1, 2), dtype=float))
        self.setattr_dataset('global_optima')

        # create data structures to store all sweep results: (v_opt_ax_0, v_opt_ax_1, v_fit_err, slope_ampl_mag)
        self.set_dataset('sweep_results', zeros((self.iterations * 2, 2, 4), dtype=float))
        self.setattr_dataset('sweep_results')

        # store most recent predicted values used
        # note: use center of voltage scan range for starting voltages
        self._host_voltage_optima_current = array([mean(self.dc_scan_range_volts_axis_0),
                                                      mean(self.dc_scan_range_volts_axis_1)])
        self._host_att_db_list = array([self.att_mode_0_db, self.att_mode_1_db], dtype=int32)

        # create data structures to hold demodulated counts
        self._host_demod_holder =       zeros((2, 100, 4), dtype=float)
        self._host_demod_holder_idx =   array([0, 0])

        # guess initial mode vectors
        # todo: improve & generalize
        if self.freq_mode_0_ftw < self.freq_mode_1_ftw:
            self.INITIAL_MODE_VECTORS = array([-1., 1.]) # [RF, EGGS]
        else:
            self.INITIAL_MODE_VECTORS = array([1., -1.]) # [EGGS, RF]

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
        # get DDS CPLD att values so ARTIQ remembers them
        self.dds_parametric.cpld.get_att_mu()

        # prepare cooling beams
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
        self.dds_parametric.set_phase_absolute()
        delay_mu(10000)
        self.dds_parametric.set_profile(self.profile_dds_parametric)
        delay_mu(8000)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        """
        Main sequence of experiment.
        """
        # predict initial optimum
        self.predict_optimum()

        # run given number of iterations
        for _iter_num in range(self.iterations):
            '''
            SWEEP VOLTAGE AXIS 0
            '''
            # prepare adapted voltage vectors and attenuations for scan
            voltage_scan_axis0_v_arr = self.prepare_voltage_scan(0)
            att_mu_mode_list = self.prepare_att(voltage_scan_axis0_v_arr)
            # scan!
            self.scan_voltage_axis(att_mu_mode_list,
                                   self.dc_channel_axis_0_num,
                                   voltage_scan_axis0_v_arr)

            # process optima to check convergence & store results
            self.process_optimum(0)
            convergence_status = self.predict_optimum()
            if convergence_status is True:
                self.core.break_realtime()
                return
            # support graceful termination
            self.check_termination()

            '''
            SWEEP VOLTAGE AXIS 1
            '''
            # prepare adapted voltage vectors and attenuations for scan
            voltage_scan_axis1_v_arr = self.prepare_voltage_scan(1)
            att_mu_mode_list = self.prepare_att(voltage_scan_axis1_v_arr)
            # scan!
            self.scan_voltage_axis(att_mu_mode_list,
                                   self.dc_channel_axis_1_num,
                                   voltage_scan_axis1_v_arr)

            # process optima to check convergence & store results
            self.process_optimum(1)
            convergence_status = self.predict_optimum()
            if convergence_status is True:
                self.core.break_realtime()
                return
            # support graceful termination
            self.check_termination()


    '''
    HELPER FUNCTIONS
    '''
    @rpc
    def predict_optimum(self) -> TBool:
        """
        Predict the location of the global micromotion optimum using
        linear fits to the measured optima of each mode.
        Returns:
            TBool:  whether we have converged on the global optimum.
        """
        print("\tBEGIN SWEEP #{:d}: {:}".format(self._host_sweep_counter,
                                                self.dc_channel_axes_names[self._host_sweep_counter % 2]))
        fitter = fitLineLinear()

        # if more than 2 minima, fit line to extract optimum
        if self._host_sweep_counter >= 2:
            # only use relatively recent results to reduce errors
            idx_min = max(self._host_sweep_counter - 4, 0)
            idx_max = self._host_sweep_counter
            optima_tmp = self.sweep_results[idx_min:idx_max, :, 0:2]

            # fit a line to the optimum for each mode
            fit_mode_0 = fitter.fit(optima_tmp[:, 0])
            fit_mode_1 = fitter.fit(optima_tmp[:, 1])

            # calculate optima as intersection of the fit lines
            opt_v_axis_0 = - (fit_mode_1[0] - fit_mode_0[0]) / (fit_mode_1[1] - fit_mode_0[1])
            opt_v_axis_1 = (fit_mode_1[1] * fit_mode_0[0] - fit_mode_0[1] * fit_mode_1[0]) / (fit_mode_1[1] - fit_mode_0[1])

        # if we only have 1 sweep, set optima as mean of previous sweep and use
        # initial vector guess to predict optima of other mode
        elif self._host_sweep_counter == 1:
            # assume voltage axis 0 was scanned
            m0, m1 = self.INITIAL_MODE_VECTORS
            b0, b1 = self.sweep_results[0, :, 0] - self.INITIAL_MODE_VECTORS * mean(self.sweep_results[0, :, 1])
            opt_v_axis_0 = (m1 * b0 - m0 * b1) / (m1 - m0)
            opt_v_axis_1 = - (b1 - b0) / (m1 - m0)

            # # old: if we only have 1 sweep, guess optima as mean of previous sweep
            # # note: doesn't propagate change to mode 1 b/c it wasn't scanned
            # opt_v_axis_0, opt_v_axis_1 = mean(self.sweep_results[0, :, 0:2], axis=0)
            # print("\n\t\tDEBUG: opt old - {} {}".format(opt_v_axis_0, opt_v_axis_1))

        # if we have no data, simply use the current voltages
        else:
            # # old: if we have no data, set optima as center of scan ranges
            # opt_v_axis_0 = mean(self.dc_scan_range_volts_list[0])
            # opt_v_axis_1 = mean(self.dc_scan_range_volts_list[1])
            opt_v_axis_0 = self.trap_dc.trap_dc.voltage(self.dc_channel_axis_0_num)
            opt_v_axis_1 = self.trap_dc.trap_dc.voltage(self.dc_channel_axis_1_num)
        print('\t\tPredicted Optimum: ({:.2f} V, {:.2f} V)'.format(opt_v_axis_0, opt_v_axis_1))

        # ensure voltages are within bounds before we set them
        if (opt_v_axis_0 < self.dc_scan_range_volts_list[0, 0]) or (opt_v_axis_0 > self.dc_scan_range_volts_list[0, 1])\
                or (opt_v_axis_1 < self.dc_scan_range_volts_list[1, 0]) or (opt_v_axis_1 > self.dc_scan_range_volts_list[1, 1]):
            raise ValueError("Error: predicted global optimum outside valid scan range.")
        # set voltages to optimum
        self.trap_dc.trap_dc.voltage_fast(self.dc_channel_axis_0_num, opt_v_axis_0)
        self.trap_dc.trap_dc.voltage_fast(self.dc_channel_axis_1_num, opt_v_axis_1)

        # store results
        opt_v_arr = array([opt_v_axis_0, opt_v_axis_1])
        self._host_voltage_optima_current = opt_v_arr
        self.mutate_dataset('global_optima', self._host_sweep_counter, opt_v_arr)

        # check successive optima for convergence
        if self._host_sweep_counter >= 2:
            if norm(opt_v_arr - self.global_optima[self._host_sweep_counter - 1]) < self.CONVERGENCE_VOLTAGE_V:
                print('STOPPING EARLY - CONVERGED ON VOLTAGE OPTIMUM (distance < {:.2f}).'.format(self.CONVERGENCE_VOLTAGE_V))
                return True

        return False

    @rpc
    def prepare_voltage_scan(self, voltage_axis: TInt32) -> TArray(TFloat, 1):
        """
        Use existing data to generate the voltage list for the sweep.
        :param voltage_axis: the voltage axis to be swept.
        :return: an array of voltages to scan.
        """
        # if we have data, use relative difference between previous optima to set the scan range
        if self._host_sweep_counter >= 1:
            # get previous optima
            opt_prev_mode0 = self.sweep_results[self._host_sweep_counter - 1, 0, 0:2]
            opt_prev_mode1 = self.sweep_results[self._host_sweep_counter - 1, 1, 0:2]
            # set voltage range depending on difference between previous optima
            # idea is to scan wide enough to hit both optima
            if voltage_axis == 0:
                voltage_range = 1.4 * abs(opt_prev_mode0[1] - opt_prev_mode1[1])
            else:
                voltage_range = 1.4 * abs(opt_prev_mode0[0] - opt_prev_mode1[0])

            # ensure voltage range is greater than 2V to get sufficient excitation
            voltage_range = max(voltage_range, 2.)

            # set voltage scan range
            voltage_max_v = self._host_voltage_optima_current[voltage_axis] + 0.5 * voltage_range
            voltage_min_v = self._host_voltage_optima_current[voltage_axis] - 0.5 * voltage_range

        # otherwise, scan about existing center, w/~40% of range about existing center
        # # old: otherwise, simply set the voltage range as half the bounds
        else:
            voltage_min_v, voltage_max_v = self.dc_scan_range_volts_list[voltage_axis]
            # voltage_center = 0.5 * (voltage_max_v + voltage_min_v)
            voltage_center = self.global_optima[self._host_sweep_counter, 0]
            voltage_range = 0.5 * (voltage_max_v - voltage_min_v)
            voltage_max_v = voltage_center + 0.5 * voltage_range
            voltage_min_v = voltage_center - 0.5 * voltage_range

        # ensure voltages are within bounds
        voltage_max_v = min(self.dc_scan_range_volts_list[voltage_axis, 1], voltage_max_v)
        voltage_min_v = max(self.dc_scan_range_volts_list[voltage_axis, 0], voltage_min_v)
        print('\t\tVoltage Scan Range (Axis {:d}, {:}): [{:.2f}, {:.2f}]'.format(voltage_axis,
                                                                                 self.dc_channel_axes_names[voltage_axis],
                                                                                 voltage_min_v, voltage_max_v))

        # create voltage scan array
        voltage_scan_arr = linspace(voltage_min_v, voltage_max_v, self.num_steps)
        shuffle(voltage_scan_arr)
        return voltage_scan_arr

    @rpc
    def prepare_att(self, voltage_scan_v_arr: TArray(TFloat, 1)) -> TList(TInt32):
        """
        Adjust attenuations to maximize the received signal as we progressively reduce the scan range.
        This function aims to keep the correlated amplitude at ~0.2 (i.e. 20% correlation).
        :param voltage_scan_v_arr: the list of voltages to scan.
        :return: the ideal attenuations to use for each mode.
        """
        # if we have data, use complex fit parameters to guess the optimal attenuation
        if self._host_sweep_counter >= 1:
            # get range of voltage scan to calculate ideal complex slope magnitude
            voltage_range = 0.5 * (max(voltage_scan_v_arr) - min(voltage_scan_v_arr))
            slope_mag_ideal = self.OPT_CORR_AMPL_FRAC / norm((voltage_range, self.GUESS_CORR_AMPL_GAMMA))

            # get most recently fitted complex slope values
            slope_mag_prev = self.sweep_results[self._host_sweep_counter - 1, :, 3]

            # calculate ideal attenuations
            att_db_change = (slope_mag_ideal - slope_mag_prev) / self.CORR_AMPL_ATT_SLOPE
            att_list_db = list(self._host_att_db_list + att_db_change)

        # otherwise, simply set the attenuations to the input values
        else:
            att_list_db = [self.att_mode_0_db, self.att_mode_1_db]

        # ensure attenuations are within bounds
        att_list_db = [min(self.att_db_bounds_list[1], att_db) for att_db in att_list_db]
        att_list_db = [max(self.att_db_bounds_list[0], att_db) for att_db in att_list_db]
        self._host_att_db_list = array(att_list_db)
        print('\t\tAttenuations:')
        print('\t\t\tMode 0 ({:.1f} kHz): {:.1f} dB\n\t\t\tMode 1 ({:.1f} kHz): {:.1f} dB'.format(
            self.freq_mode_0_khz, att_list_db[0],
            self.freq_mode_1_khz, att_list_db[1])
        )

        # convert results to machine units and return
        att_list_mu = [self.dds_parametric.cpld.att_to_mu(att_db * dB) for att_db in att_list_db]
        return att_list_mu

    @kernel(flags={"fast-math"})
    def scan_voltage_axis(self, mode_att_mu_list: TList(TInt32), dc_channel_num: TInt32,
                          voltage_scan_v_arr: TArray(TFloat, 1)) -> TNone:
        """
        Characterize micromotion along a voltage axis via parametric excitation.
        :param mode_att_mu_list: list of attenuations to use (in machine units).
        :param dc_channel_num: the voltage channel number to be scanned.
        :param voltage_scan_v_arr: the list of voltages to scan.
        """
        # scan voltage configurations in the voltage vector
        for voltage_v in voltage_scan_v_arr:

            # set DC voltage and synchronize hardware clock with timeline
            self.trap_dc.trap_dc.voltage_fast(dc_channel_num, voltage_v)
            self.core.break_realtime() # add slack so now_mu() isn't in the past
            self.core.wait_until_mu(now_mu())
            # add extra delay for voltages to settle
            delay_mu(self.time_dc_synchronize_delay_mu)

            # run N reps at target voltage (since voltage adjustment dominates overheads)
            for rep_num in range(self.repetitions_per_voltage):
                # swap between modes within each rep
                for mode_idx in self.mode_list_idx:

                    # prepare hardware for parametric demodulation
                    mode_freq_ftw = self.freq_mode_ftw_list[mode_idx]
                    mode_att_mu =   mode_att_mu_list[mode_idx]

                    self.core.break_realtime()
                    self.dds_parametric.set_att_mu(mode_att_mu)
                    self.dds_parametric.set_mu(mode_freq_ftw, asf=self.dds_parametric.ampl_modulation_asf,
                                               profile=self.profile_dds_parametric,
                                               phase_mode=PHASE_MODE_CONTINUOUS)
                    self.pump.set_profile(self.profile_397_parametric)
                    delay_mu(25000)

                    # run parametric excitation and get timestamps
                    pmt_timestamp_list = self.parametric_subsequence.run()

                    # clean up and process results (e.g. demodulate counts)
                    self.pump.rescue()  # leave beams on rescue while idle
                    self.pmt.clear_inputs()
                    self._demodulate_counts(mode_idx, mode_freq_ftw, voltage_v, pmt_timestamp_list)
                    self.check_termination()

    @rpc
    def _demodulate_counts(self, mode_idx: TInt32, freq_mu: TInt32, voltage_v: TFloat,
                          timestamp_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Demodulate the PMT timestamps with respect to the parametric modulation frequency.
        :param mode_idx: the voltage axis being swept.
        :param freq_ftw: the modulation frequency (as a 32-bit frequency tuning word).
        :param voltage_v: the current shim voltage (in volts).
        :param timestamp_mu_list: the list of timestamps (in machine units) to demodulate.
        """
        # convert arguments to human units
        freq_khz =      self.dds_parametric.ftw_to_frequency(freq_mu) / kHz
        timestamps_s =  self.core.mu_to_seconds(array(timestamp_mu_list))
        
        # digitally demodulate counts and convert correlated signal to polar coordinates
        correlated_signal = mean(exp((2.j * pi * freq_khz * 1e3) * timestamps_s))
        correlated_ampl =   abs(correlated_signal)
        correlated_phase =  angle(correlated_signal)
        # extract count rate in seconds
        count_rate_hz =     len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # store results in global parent dataset for post-processing
        res_arr = array([freq_khz, voltage_v, correlated_ampl, correlated_phase])
        self.update_results(*res_arr, count_rate_hz)
        # store results in holder for immediate processing
        self._host_demod_holder[mode_idx, self._host_demod_holder_idx[mode_idx]] = res_arr
        self._host_demod_holder_idx[mode_idx] += 1

    @rpc
    def process_optimum(self, voltage_axis: TInt32) -> TNone:
        """
        Process demodulated counts to extract voltage optima for each mode.
        Necessary for kernel to interface with results since raw data is store on host-side.
        Note: this should only be called at the end of a voltage scan, since this updates internal
            counters/iterators (i.e. this function defines the end of a voltage scan).
        # todo: check results to ensure max corr ampl is above a given threshold
        :param voltage_axis: the voltage axis that was swept.
        """
        opt_res_store = []

        # process optima for each mode
        for mode_idx in self.mode_list_idx:

            # get all recent mode data
            res_arr = self._host_demod_holder[mode_idx, :self._host_demod_holder_idx[mode_idx]]

            # sort & group array by voltage (column 1), then average demodulated signal
            # note: sort is necessary here since groupby requires key/values to be contiguous
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
            print("\t\tMode {:d} Opt. (Axis {:d}, {:}): {:.2f} +/- {:.3f} V".format(
                mode_idx, voltage_axis, self.dc_channel_axes_names[voltage_axis],
                opt_voltage_v, opt_voltage_err)
            )

            # check optima for errors
            if ((opt_voltage_v < self.dc_scan_range_volts_list[voltage_axis, 0]) or
                    (opt_voltage_v > self.dc_scan_range_volts_list[voltage_axis, 1])):
                raise ValueError("Error: Mode {:d} voltage out of range: {:f} V.".format(
                    mode_idx, opt_voltage_v))
            elif abs(opt_voltage_err) > 2.0:
                raise ValueError("Error: Mode {:d} optimum uncertainty exceeds bounds: {:f} V.".format(
                    mode_idx, opt_voltage_err))

            # get full optima vector
            if voltage_axis == 0:
                opt_v_mode1 = self._host_voltage_optima_current[1]
                opt_v_vector = array([opt_voltage_v, opt_v_mode1])
            elif voltage_axis == 1:
                opt_v_mode0 = self._host_voltage_optima_current[0]
                opt_v_vector = array([opt_v_mode0, opt_voltage_v])

            # store results
            opt_res_store.append(array([*opt_v_vector, opt_voltage_err, slope_ampl_mag]))

        # store results in dataset
        self.mutate_dataset('sweep_results', self._host_sweep_counter, opt_res_store)

        # clear results and update indices
        self._host_demod_holder[:] = 0.
        self._host_demod_holder_idx[:] = 0
        self._host_sweep_counter += 1


    '''
    ANALYZE EXPERIMENT
    '''
    def analyze_experiment(self):
        # todo: analyze
        pass

