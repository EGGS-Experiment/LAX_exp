import numpy as np
from itertools import groupby

from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.system.subsequences import ParametricExcite, RescueIon
import LAX_exp.experiments.parametric.ParametricSweep as ParametricSweep

import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config

# todo: set kernel invariants


class MicromotionCompensationRDX(ParametricSweep.ParametricSweep, Experiment):
    """
    Experiment: Micromotion Compensation RDX

    Characterize the micromotion along a mode by applying a parametric excitation on the
    trap RF while scanning shim voltages, then attempt to algorithmically compensate for it.
    """
    name = 'Micromotion Compensation RDX'

    kernel_invariants = {
        "freq_mod0_ftw", "freq_mod1_ftw", "freq_mode_ftw_list", "mode_list_idx",
        "att_mod0_mu", "att_mod1_mu",
        "dc_channel_axis_0_num", "dc_channel_axis_1_num", "dc_scan_range_volts_list", "time_dc_synchronize_delay_mu",
        "ampl_cooling_asf", "freq_cooling_ftw", "time_cooling_holdoff_mu",
        "cxn", "dc"
    }


    def build_experiment(self):
        # get DC channel configuration dictionary
        self.dc_config_channeldict =                                dc_config.channeldict

        # core arguments
        self.setattr_argument("iterations",                 NumberValue(default=5, ndecimals=0, step=1, min=1, max=10))

        # general configuration
        self.setattr_argument("repetitions_per_voltage",    NumberValue(default=2, ndecimals=0, step=1, min=1, max=100), group='configuration')
        self.setattr_argument("num_steps",                  NumberValue(default=10, ndecimals=0, step=1, min=5, max=100), group='configuration')
        self.setattr_argument("adaptive",                   BooleanValue(default=True), group='configuration')

        # modulation - mode #1
        self.setattr_argument("freq_mod0_khz",              NumberValue(default=1274.13, ndecimals=3, step=10, min=1, max=10000), group='modulation')
        self.setattr_argument("att_mod0_db",                NumberValue(default=20, ndecimals=1, step=0.5, min=0, max=31.5), group='modulation')
        # modulation - mode #2
        self.setattr_argument("freq_mod1_khz",              NumberValue(default=1565.92, ndecimals=3, step=10, min=1, max=10000), group='modulation')
        self.setattr_argument("att_mod1_db",                NumberValue(default=15, ndecimals=1, step=0.5, min=0, max=31.5), group='modulation')

        # shim voltages
        self.setattr_argument("dc_channel_axis_0",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='V Shim'), group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_0", PYONValue([45, 85]), group='voltages')
        self.setattr_argument("dc_channel_axis_1",          EnumerationValue(list(self.dc_config_channeldict.keys()), default='H Shim'), group='voltages')
        self.setattr_argument("dc_scan_range_volts_axis_1", PYONValue([30, 70]), group='voltages')

        # cooling
        self.setattr_argument("ampl_cooling_pct",           NumberValue(default=24, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",           NumberValue(default=105, ndecimals=6, step=1, min=1, max=500), group='cooling')

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_parametric')

        # get relevant subsequences
        self.parametric_subsequence =                           ParametricExcite(self)
        self.rescue_subsequence =                               RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare the experiment.
        """
        '''
        VOLTAGES
        '''
        # sanitize voltage input
        for voltage_arr in [self.dc_scan_range_volts_axis_0, self.dc_scan_range_volts_axis_1]:
            # check voltage scan ranges have correct format
            if (type(voltage_arr) is not list) or (len(voltage_arr) != 2):
                raise Exception("InputError: voltage scan ranges have incorrect type.")

            # check voltages are all in valid range
            if any([(voltage < 0) or (voltage > 110) for voltage in voltage_arr]):
                raise Exception("InputError: voltage range is out of bounds.")

        # get DC channel numbers
        self.dc_channel_axis_0_num =    self.dc_config_channeldict[self.dc_channel_axis_0]['num']
        self.dc_channel_axis_1_num =    self.dc_config_channeldict[self.dc_channel_axis_1]['num']

        # store voltage bounds conveniently for programmatic access
        self.dc_scan_range_volts_list = np.array([np.sort(self.dc_scan_range_volts_axis_0),
                                                  np.sort(self.dc_scan_range_volts_axis_1)])

        # holdoff period after we set a voltage to allow it to settle
        self.time_dc_synchronize_delay_mu = self.core.seconds_to_mu(988 * ms)


        '''
        MODULATION
        '''
        # convert modulation (mode #0) parameters to machine units
        self.freq_mod0_ftw =            self.dds_parametric.frequency_to_ftw(self.freq_mod0_khz * kHz)
        self.att_mod0_mu =              att_to_mu(self.att_mod0_db * dB)

        # convert modulation (mode #1) parameters to machine units
        self.freq_mod1_ftw =            self.dds_parametric.frequency_to_ftw(self.freq_mod1_khz * kHz)
        self.att_mod1_mu =              att_to_mu(self.att_mod1_db * dB)

        # store modulation values conveniently for programmatic access
        self.freq_mode_ftw_list =       np.array([self.freq_mod0_ftw, self.freq_mod1_ftw], dtype=np.int32)
        self.mode_list_idx =            list(range(len(self.freq_mode_ftw_list)))


        '''
        COOLING
        '''
        # convert cooling parameters to machine units
        self.ampl_cooling_asf =                 self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.freq_cooling_ftw =                 self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.time_cooling_holdoff_mu =          self.core.seconds_to_mu(3 * ms)


        '''
        DATA STRUCTURES
        '''
        # tmp remove
        self.set_dataset('global_optima', np.zeros((self.iterations * 2, 2), dtype=float))
        self.setattr_dataset('global_optima')
        # tmp remove

        # create data structures to store voltage optima
        self.set_dataset('dc_voltage_optima', np.zeros((self.iterations * 2, 2, 3), dtype=float))
        self.setattr_dataset('dc_voltage_optima')

        self._host_opt_holder =         np.zeros((self.iterations * 2, 2, 2), dtype=float)
        self._host_sweep_counter =      0

        # note: use center of voltage scan range for starting voltages
        self._host_voltage_optima_current = np.array([np.mean(self.dc_scan_range_volts_axis_0),
                                                      np.mean(self.dc_scan_range_volts_axis_1)])

        # create data structures to hold demodulated counts
        self._host_demod_holder =       np.zeros((2, 100, 4), dtype=float)
        self._host_demod_holder_idx =   np.array([0, 0])


        '''
        LABRAD
        '''
        # connect to labrad and get DC server
        self.cxn =  labrad.connect(environ['LABRADHOST'],
                                   port=7682, tls_mode='off',
                                   username='', password='lab')
        self.dc =   self.cxn.dc_server

    @property
    def results_shape(self):
        # x2 for 2 voltage axes; x2 for 2 modes; x2 for safety
        return (self.iterations * 2 * self.num_steps * self.repetitions_per_voltage * 2 * 2,
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # set up labrad devices via RPC
        self.prepareDevicesLabrad()
        self.core.break_realtime()
        
        # get DDS CPLD att values so ARTIQ remembers them
        self.dds_parametric.cpld.get_att_mu()
        self.core.break_realtime()

        # set cooling beams
        self.pump.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.pump.set_profile(0)
        self.pump.on()
        self.repump_cooling.on()
        self.repump_qubit.on()

        # set up DDS for modulation
        self.dds_parametric.set_phase_absolute()


    @kernel(flags={"fast-math"})
    def run_main(self):
        """
        Main sequence of experiment.
        """
        # run given number of iterations
        for _iter_num in range(self.iterations):
            self.core.break_realtime()

            '''
            PREPARE OPTIMIZATION LOOP
            '''
            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


            '''
            SWEEP VOLTAGE AXIS 0
            '''
            # calculate voltage and set optima
            self.predict_optimum()
            self.core.break_realtime()

            # generate voltage vector array and scan
            voltage_scan_axis0_v_arr = self.prepare_voltage_scan(0)
            self.core.break_realtime()
            self.scan_voltage_axis([self.att_mod0_mu, self.att_mod1_mu],
                                   self.dc_channel_axis_0_num,
                                   voltage_scan_axis0_v_arr)
            self.core.break_realtime()

            # process and store results
            self.process_optimum(0)
            self.core.break_realtime()


            '''
            SWEEP VOLTAGE AXIS 1
            '''
            # calculate voltage and set optima
            self.predict_optimum()
            self.core.break_realtime()

            # generate voltage vector array and scan
            voltage_scan_axis1_v_arr = self.prepare_voltage_scan(1)
            self.core.break_realtime()
            self.scan_voltage_axis([self.att_mod0_mu, self.att_mod1_mu],
                                   self.dc_channel_axis_1_num,
                                   voltage_scan_axis1_v_arr)
            self.core.break_realtime()

            # process and store results
            self.process_optimum(1)
            self.core.break_realtime()



    '''
    HELPER FUNCTIONS
    '''
    @rpc
    def predict_optimum(self) -> TNone:
        """
        Predict the location of the global micromotion optimum using
        linear fits to the measured optima of each mode.
        """
        # tmp remove
        print("\t\tbegin predict optima")
        # tmp remove

        # if more than 2 minima, fit line to extract optimum
        if self._host_sweep_counter >= 2:
            # only use relatively recent results to reduce errors
            idx_min = max(self._host_sweep_counter - 4, 0)
            idx_max = self._host_sweep_counter
            optima_tmp = self._host_opt_holder[idx_min:idx_max, :, :]

            # tmp remove
            print(optima_tmp)
            # tmp remove

            # fit a line to the optimum for each mode
            fit_mode_0 = fitLineLinear(optima_tmp[:, 0])
            fit_mode_1 = fitLineLinear(optima_tmp[:, 1])

            # calculate optima as intersection of the fit lines
            opt_v_axis_0 = (fit_mode_0[1] * fit_mode_1[0] + fit_mode_0[0]) / (1. - fit_mode_0[1] * fit_mode_1[1])
            opt_v_axis_1 = (fit_mode_1[1] * fit_mode_0[0] + fit_mode_1[0]) / (1. - fit_mode_0[1] * fit_mode_1[1])

        # if we only have 1 sweep, set optima as mean of previous sweep
        elif self._host_sweep_counter == 1:
            opt_v_axis_0, opt_v_axis_1 = np.mean(self._host_opt_holder[0, :, :], axis=0)

        # if we have no data, then set optima as center of scan ranges
        else:
            opt_v_axis_0 = np.mean(self.dc_scan_range_volts_list[0])
            opt_v_axis_1 = np.mean(self.dc_scan_range_volts_list[1])
        print('\t\tPredicted Optimum: [{:.2f} V, {:.2f} V]'.format(opt_v_axis_0, opt_v_axis_1))

        # ensure voltages are within bounds
        if (opt_v_axis_0 < self.dc_scan_range_volts_list[0, 0]) or (opt_v_axis_0 > self.dc_scan_range_volts_list[0, 1])\
                or (opt_v_axis_1 < self.dc_scan_range_volts_list[1, 0]) or (opt_v_axis_1 > self.dc_scan_range_volts_list[1, 1]):
            raise Exception("Error: global optimum predicted to be outside valid scan range.")

        # set voltages to optimum
        self.voltage_set(self.dc_channel_axis_0_num, opt_v_axis_0)
        self.voltage_set(self.dc_channel_axis_1_num, opt_v_axis_1)

        # store results
        opt_v_arr = np.array([opt_v_axis_0, opt_v_axis_1])
        self._host_voltage_optima_current = opt_v_arr
        self.mutate_dataset('global_optima', self._host_sweep_counter, opt_v_arr)

        # tmp remove
        print("\t\tend predict optima")
        # tmp remove

    @rpc
    def prepare_voltage_scan(self, voltage_axis: TInt32) -> TArray(TFloat, 1):
        """
        # todo: document
        Arguments:
            voltage_axis    (TInt32)    : the voltage axis to be swept.
        Returns:
            TArray(TFloat, 1)   : an array of voltages to scan.
        """
        print("\tBEGIN SWEEP: {:d}".format(self._host_sweep_counter))

        # if we have data, use relative difference between previous optima to set the scan range
        if self._host_sweep_counter >= 1:
            # get previous optima
            opt_prev_mode0 = self._host_opt_holder[self._host_sweep_counter - 1, 0]
            opt_prev_mode1 = self._host_opt_holder[self._host_sweep_counter - 1, 1]
            if voltage_axis == 0:
                voltage_range = 1.4 * abs(opt_prev_mode0[1] - opt_prev_mode1[1])
            else:
                voltage_range = 1.4 * abs(opt_prev_mode0[0] - opt_prev_mode1[0])

            # ensure voltage range is greater than 2V to get sufficient excitation
            voltage_range = max(voltage_range, 2.)

            # set voltage scan range
            voltage_max_v = self._host_voltage_optima_current[voltage_axis] + 0.5 * voltage_range
            voltage_min_v = self._host_voltage_optima_current[voltage_axis] - 0.5 * voltage_range

        # otherwise, simply set the voltage range as half the bounds
        else:
            voltage_min_v, voltage_max_v = self.dc_scan_range_volts_list[voltage_axis]
            voltage_center = 0.5 * (voltage_max_v + voltage_min_v)
            voltage_range = 0.5 * (voltage_max_v - voltage_min_v)
            voltage_max_v = voltage_center + 0.5 * voltage_range
            voltage_min_v = voltage_center - 0.5 * voltage_range

        # todo: adjust attenuation based on correlated amplitude

        # ensure voltages are within bounds
        voltage_max_v = min(self.dc_scan_range_volts_list[voltage_axis, 1], voltage_max_v)
        voltage_min_v = max(self.dc_scan_range_volts_list[voltage_axis, 0], voltage_min_v)
        # todo: print name IN ADDITION TO axis
        print('\t\tVoltage Scan Range (Axis {:d}): [{:.2f}, {:.2f}]'.format(voltage_axis, voltage_min_v, voltage_max_v))

        # create voltage scan array
        voltage_scan_arr = np.linspace(voltage_min_v, voltage_max_v, self.num_steps)
        np.random.shuffle(voltage_scan_arr)
        return voltage_scan_arr

    @rpc
    def prepare_att(self):
        pass

    @kernel(flags={"fast-math"})
    def scan_voltage_axis(self, mode_att_mu_list: TList(TInt32), dc_channel_num: TInt32,
                          voltage_scan_v_arr: TArray(TFloat, 1)) -> TNone:
        """
        Characterize micromotion along a voltage axis via parametric excitation.
        Arguments:
            mode_att_mu_list    (TList(TInt32))     : list of attenuations to use (in machine units).
            dc_channel_num      (TInt32)            : the voltage channel number to be scanned.
            voltage_scan_v_arr  (TArray(TFloat, 1)) : the list of voltages to scan.
        """
        # iterate over repetitions per voltage
        for rep_num in range(self.repetitions_per_voltage):

            # scan voltage configurations in the voltage vector
            for voltage_v in voltage_scan_v_arr:

                # set DC voltage and synchronize hardware clock with timeline
                self.voltage_set(dc_channel_num, voltage_v)
                self.core.wait_until_mu(now_mu())
                # add extra delay for voltages to settle
                delay_mu(self.time_dc_synchronize_delay_mu)

                # get parametric excitation data for both modes
                for mode_idx in self.mode_list_idx:

                    # extract values
                    mode_freq_ftw = self.freq_mode_ftw_list[mode_idx]
                    mode_att_mu =   mode_att_mu_list[mode_idx]

                    # prepare modulation DDS
                    self.dds_parametric.set_att_mu(mode_att_mu)
                    self.dds_parametric.set_mu(mode_freq_ftw, asf=self.dds_parametric.ampl_modulation_asf,
                                               profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
                    self.core.break_realtime()

                    # run parametric excitation and get timestamps
                    pmt_timestamp_list = self.parametric_subsequence.run()

                    # demodulate counts and store
                    self._demodulate_counts(mode_idx, mode_freq_ftw, voltage_v, pmt_timestamp_list)
                    self.core.reset()

            # rescue ion as needed
            self.rescue_subsequence.run(rep_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

    @rpc
    def _demodulate_counts(self, mode_idx: TInt32, freq_mu: TInt32, voltage_v: TFloat,
                          timestamp_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Demodulate the PMT timestamps with respect to the parametric modulation frequency.
        Arguments:
            mode_idx            (TInt32)            : the voltage axis being swept.
            freq_ftw            (TInt32)            : the modulation frequency (as a 32-bit frequency tuning word).
            voltage_v           (TFloat)            : the current shim voltage (in volts).
            timestamp_mu_list   (TArray(TInt64, 1)) : the list of timestamps (in machine units) to demodulate.
        """
        # convert arguments to human units
        freq_khz =      self.dds_parametric.ftw_to_frequency(freq_mu) / kHz
        timestamps_s =  self.core.mu_to_seconds(np.array(timestamp_mu_list))
        
        # digitally demodulate counts and convert correlated signal to polar coordinates
        correlated_signal = np.mean(np.exp((2.j * np.pi * freq_khz * 1e3) * timestamps_s))
        correlated_ampl =   np.abs(correlated_signal)
        correlated_phase =  np.angle(correlated_signal)
        # extract count rate in seconds
        count_rate_hz =     len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # store results in global parent dataset for post-processing
        res_arr = np.array([freq_khz, voltage_v, correlated_ampl, correlated_phase])
        self.update_results(*res_arr, count_rate_hz)
        # store results in holder for immediate processing
        self._host_demod_holder[mode_idx, self._host_demod_holder_idx[mode_idx]] = res_arr
        self._host_demod_holder_idx[mode_idx] += 1

    @rpc
    def process_optimum(self, voltage_axis: TInt32) -> TNone:
        """
        Retrieve recent experiment data
        Necessary for kernel to interface with results since raw data is store on host-side.
        Arguments:
            voltage_axis    (TInt32): the voltage axis that was swept.
        """
        # tmp remove
        print("\t\tbegin process optimum")
        # tmp remove

        opt_res_store = []

        # process optima for each mode
        for mode_idx in self.mode_list_idx:

            # get all recent mode data
            res_arr = self._host_demod_holder[mode_idx, :self._host_demod_holder_idx[mode_idx]]

            # sort & group array by voltage (column 1), then average demodulated signal
            # note: sort is necessary here since groupby requires key/values to be contiguous
            res_arr_tmp = res_arr[np.argsort(res_arr[:, 1])]
            res_arr_tmp = np.array([
                np.mean(np.array([val for val in group]), axis=0)
                for key, group in groupby(res_arr_tmp, lambda arr: arr[1])
            ])

            # format results into a 2D array with complex type for complex linear fitting
            results_tmp = np.array([
                res_arr_tmp[:, 1],
                res_arr_tmp[:, 2] * np.exp(1.j * res_arr_tmp[:, 3])
            ], dtype='complex128').transpose()
            # sort results by voltage (again)
            results_tmp = results_tmp[np.argsort(results_tmp[:, 1], axis=0)]

            # extract minimum mode voltage
            opt_voltage_v, opt_voltage_err = complexLinearFitMinimize(results_tmp)

            # check optima for errors
            if ((opt_voltage_v < self.dc_scan_range_volts_list[voltage_axis, 0]) or
                    (opt_voltage_v > self.dc_scan_range_volts_list[voltage_axis, 1])):
                raise Exception("Error: Mode {:d} voltage out of range: {:f} V.".format(mode_idx, opt_voltage_v))
            elif abs(opt_voltage_err) > 2.0:
                raise Exception("Error: Mode {:d} optimum uncertainty exceeds bounds: {:f} V.".format(mode_idx, opt_voltage_err))

            # get full optima vector
            if voltage_axis == 0:
                opt_v_mode1 = self._host_voltage_optima_current[1]
                opt_v_vector = np.array([opt_voltage_v, opt_v_mode1])
            elif voltage_axis == 1:
                opt_v_mode0 = self._host_voltage_optima_current[0]
                opt_v_vector = np.array([opt_v_mode0, opt_voltage_v])

            # store results in host holders
            # todo: say which voltage rod is being used instead of e.g. axis0
            print("\t\tMode {:d} Opt. (Axis {:d}): {:.2f} +/- {:.3f} V".format(mode_idx, voltage_axis, opt_voltage_v, opt_voltage_err))
            self._host_opt_holder[self._host_sweep_counter, mode_idx, :] = opt_v_vector
            opt_res_store.append(np.array([*opt_v_vector, opt_voltage_err]))

        # store results in dataset
        self.mutate_dataset('dc_voltage_optima', self._host_sweep_counter, opt_res_store)

        # clear results and update indices
        self._host_demod_holder[:] = 0.
        self._host_demod_holder_idx[:] = 0
        self._host_sweep_counter += 1

        # tmp remove
        print("\t\tend process optimum")
        # tmp remove


    # ANALYZE EXPERIMENT
    def analyze_experiment(self):
        pass

