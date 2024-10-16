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


class MicromotionCompensation(ParametricSweep.ParametricSweep, Experiment):
    """
    Experiment: Micromotion Compensation

    Characterize the micromotion along a mode by applying a parametric excitation on the
    trap RF while scanning shim voltages, then attempt to algorithmically compensate for it.
    """
    name = 'Micromotion Compensation'


    def build_experiment(self):
        # get DC channel configuration dictionary
        self.dc_config_channeldict =                                dc_config.channeldict

        # core arguments
        self.setattr_argument("iterations",                         NumberValue(default=5, precision=0, step=1, min=1, max=100))

        # modulation - general
        self.setattr_argument("repetitions_per_voltage",    NumberValue(default=2, precision=0, step=1, min=1, max=100), group='modulation_general')
        self.setattr_argument("num_steps",                  NumberValue(default=10, precision=0, step=1, min=5, max=100), group='modulation_general')
        # self.setattr_argument("max_att_mod_db",                   NumberValue(default=10, precision=1, step=0.5, min=0, max=31.5), group='modulation_general')
        self.setattr_argument("adaptive",                   BooleanValue(default=True), group='modulation_general')

        # modulation - mode #1
        self.setattr_argument("freq_mod0_khz",              NumberValue(default=1394.3, precision=3, step=10, min=1, max=10000), group='modulation_0')
        self.setattr_argument("att_mod0_db",                NumberValue(default=25, precision=1, step=0.5, min=0, max=31.5), group='modulation_0')
        self.setattr_argument("dc_channel_mod0",            EnumerationValue(list(self.dc_config_channeldict.keys()), default='V Shim'), group='modulation_0')
        self.setattr_argument("dc_scan_range_volts_mod0",   PYONValue([55, 75]), group='modulation_0')

        # modulation - mode #2
        self.setattr_argument("freq_mod1_khz",              NumberValue(default=1634., precision=3, step=10, min=1, max=10000), group='modulation_1')
        self.setattr_argument("att_mod1_db",                NumberValue(default=19, precision=1, step=0.5, min=0, max=31.5), group='modulation_1')
        self.setattr_argument("dc_channel_mod1",            EnumerationValue(list(self.dc_config_channeldict.keys()), default='H Shim'), group='modulation_1')
        self.setattr_argument("dc_scan_range_volts_mod1",   PYONValue([30, 55]), group='modulation_1')

        # cooling
        self.setattr_argument("ampl_cooling_pct",           NumberValue(default=30, precision=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",           NumberValue(default=105, precision=6, step=1, min=1, max=500), group='cooling')

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
        for voltage_arr in [self.dc_scan_range_volts_mod0, self.dc_scan_range_volts_mod1]:
            # check voltage scan ranges have correct format
            if (type(voltage_arr) is not list) or (len(voltage_arr) != 2):
                raise Exception("InputError: voltage scan ranges have incorrect type.")

            # check voltages are all in valid range
            if any([(voltage < 0) or (voltage > 110) for voltage in voltage_arr]):
                raise Exception("InputError: voltage range is out of bounds.")

        # convert modulation (mode #0) parameters to machine units
        self.freq_mod0_ftw =                    self.dds_parametric.frequency_to_ftw(self.freq_mod0_khz * kHz)
        self.att_mod0_mu =                      att_to_mu(self.att_mod0_db * dB)
        self.dc_channel_mod0_num =              self.dc_config_channeldict[self.dc_channel_mod0]['num']

        # convert modulation (mode #1) parameters to machine units
        self.freq_mod1_ftw =                    self.dds_parametric.frequency_to_ftw(self.freq_mod1_khz * kHz)
        self.att_mod1_mu =                      att_to_mu(self.att_mod1_db * dB)
        self.dc_channel_mod1_num =              self.dc_config_channeldict[self.dc_channel_mod1]['num']

        # store scan range bounds
        self.dc_scan_range_volts_list =         np.array([np.sort(self.dc_scan_range_volts_mod0),
                                                          np.sort(self.dc_scan_range_volts_mod1)])

        '''
        COOLING
        '''
        self.time_dc_synchronize_delay_mu =     self.core.seconds_to_mu(488 * ms)

        # convert cooling parameters to machine units
        self.ampl_cooling_asf =                 self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.freq_cooling_ftw =                 self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.time_cooling_holdoff_mu =          self.core.seconds_to_mu(3 * ms)

        '''
        DATA STRUCTURES
        '''
        # create array to store voltage optima
        self.set_dataset('dc_voltage_optima_0', np.zeros((self.iterations, 3), dtype=float))
        self.set_dataset('dc_voltage_optima_1', np.zeros((self.iterations, 3), dtype=float))
        self.setattr_dataset('dc_voltage_optima_0')
        self.setattr_dataset('dc_voltage_optima_1')

        # create array to store current voltage values
        # note: use median of voltage scan range for starting voltages
        self.dc_voltage_optima_running =        np.array([np.mean(self.dc_scan_range_volts_mod0),
                                                          np.mean(self.dc_scan_range_volts_mod1)])
        
        # create data strutures to hold intermediate results
        self._res_holder_tmp = np.zeros((100, 3), dtype=float)
        self._res_holder_iter = 0

        '''
        LABRAD
        '''
        # connect to labrad and get DC server
        self.cxn =                              labrad.connect(environ['LABRADHOST'],
                                                               port=7682, tls_mode='off',
                                                               username='', password='lab')
        self.dc =                               self.cxn.dc_server

    @property
    def results_shape(self):
        return (self.iterations * self.repetitions_per_voltage * self.num_steps * 2,
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        with parallel:
            # set cooling beams
            with sequential:
                self.pump.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
                self.pump.set_profile(0)
                self.pump.on()
                self.repump_cooling.on()
                self.repump_qubit.on()

            # set up DDS for modulation
            self.dds_parametric.set_phase_absolute()

            # set up labrad devices via RPC
            self.prepareDevicesLabrad()
            self.core.break_realtime()


    def run_main(self):
        # run given number of iterations
        for _iter_num in range(self.iterations):

            '''
            PREPARE OPTIMIZATION LOOP
            '''
            # support graceful termination
            self.check_termination()

            # print iteration message
            print("\n\tBEGIN: ITERATION #{:d}".format(_iter_num))

            # calculate voltage and set optima
            self.dc_voltage_optima_running = self._predict_optima(_iter_num)
            self.voltage_set(self.dc_channel_mod0_num, self.dc_voltage_optima_running[0])
            self.voltage_set(self.dc_channel_mod1_num, self.dc_voltage_optima_running[1])

            '''
            OPTIMIZE MODE 0
            '''
            # generate voltage vector array (mode 0) and prepare
            voltage_scan_mod0_v_arr = self._prepare_voltage_scan(0, _iter_num)
            # scan voltage on given mode
            self._scan_voltage_mode(self.freq_mod0_ftw,
                                     self.att_mod0_mu,
                                     self.dc_channel_mod0_num,
                                     voltage_scan_mod0_v_arr)
            opt_voltage_mod0, opt_err_mod0 = self._extract_optimum(self._res_holder_tmp)

            # store result
            print("\t\tMode 0 Opt: {:.2f} +/- {:.3f}".format(opt_voltage_mod0, opt_err_mod0))
            self.mutate_dataset('dc_voltage_optima_0', _iter_num, np.array([self.dc_voltage_optima_running[1],
                                                                            opt_voltage_mod0,
                                                                            opt_err_mod0]))

            # ensure optimum voltage is in range
            if (opt_voltage_mod0 < self.dc_scan_range_volts_list[0, 0]) or (opt_voltage_mod0 > self.dc_scan_range_volts_list[0, 1]):
                print("Error: Mode 0 voltage out of range.")
                return
            # ensure optimum voltage has been extracted reasonably
            elif (abs(opt_err_mod0) > 2.0):
                print("Error: unable to extract Mode 0 optimum with reasonable certainty.")
                return

            # todo: adjust attenuation based on correlated amplitude

            # set optimum voltage
            self.voltage_set(self.dc_channel_mod0_num, opt_voltage_mod0)
            self.dc_voltage_optima_running[0] = opt_voltage_mod0


            '''
            OPTIMIZE MODE 1
            '''
            # generate voltage vector array (mode 0) and prepare
            voltage_scan_mod1_v_arr = self._prepare_voltage_scan(1, _iter_num)
            # scan voltage on given mode
            self._scan_voltage_mode(self.freq_mod1_ftw,
                                     self.att_mod1_mu,
                                     self.dc_channel_mod1_num,
                                     voltage_scan_mod1_v_arr)
            opt_voltage_mod1, opt_err_mod1 = self._extract_optimum(self._res_holder_tmp)

            # store result
            print("\t\tMode 1 Opt: {:.2f} +/- {:.3f}".format(opt_voltage_mod1, opt_err_mod1))
            self.mutate_dataset('dc_voltage_optima_1', _iter_num, np.array([self.dc_voltage_optima_running[0],
                                                                            opt_voltage_mod1,
                                                                            opt_err_mod1]))

            # ensure optimum voltage is in range
            if (opt_voltage_mod1 < self.dc_scan_range_volts_mod1[0]) or (opt_voltage_mod1 > self.dc_scan_range_volts_mod1[1]):
                print("\tError: Mode 1 voltage out of range.")
                return
            # ensure optimum voltage has been extracted reasonably
            elif (abs(opt_err_mod1) > 2.0):
                print("\tError: unable to extract Mode 1 optimum with reasonable certainty.")
                return

            # set optimum voltage
            self.voltage_set(self.dc_channel_mod1_num, opt_voltage_mod1)
            self.dc_voltage_optima_running[1] = opt_voltage_mod1


    '''
    HELPER FUNCTIONS
    '''
    @rpc
    def _predict_optima(self, iter_num: TInt32) -> TArray(TFloat, 1):
        """
        # todo: document
        Arguments:
            mode_num    (int32)         : the modulation mode.
        Returns:
                        TArray(Tfloat)  : the voltage vector of the predicted optimum.
        """
        opt0_volts, opt1_volts = (0., 0.)

        if iter_num >= 2:
            # fit a line to the optimum for each mode
            fit_mode0 = fitLineLinear(self.dc_voltage_optima_0[:iter_num, :2])
            fit_mode1 = fitLineLinear(self.dc_voltage_optima_1[:iter_num, :2])

            # calculate optima as intersection of the fit lines
            opt0_volts = (fit_mode0[1] * fit_mode1[0] + fit_mode0[0]) / (1. - fit_mode0[1] * fit_mode1[1])
            opt1_volts = (fit_mode1[1] * fit_mode0[0] + fit_mode1[0]) / (1. - fit_mode0[1] * fit_mode1[1])
        else:
            # guess the optimum as the mean of the scan range
            opt0_volts = np.mean(self.dc_scan_range_volts_list[0])
            opt1_volts = np.mean(self.dc_scan_range_volts_list[1])

        print('\tOptima: [{:.2f} V, {:.2f} V]'.format(opt0_volts, opt1_volts))


        # ensure voltages are within bounds
        opt0_volts = max(self.dc_scan_range_volts_list[0, 0], opt0_volts)
        opt0_volts = min(self.dc_scan_range_volts_list[0, 1], opt0_volts)
        opt1_volts = max(self.dc_scan_range_volts_list[1, 0], opt1_volts)
        opt1_volts = min(self.dc_scan_range_volts_list[1, 1], opt1_volts)
        return np.array([opt0_volts, opt1_volts])

    @rpc
    def _prepare_voltage_scan(self, mode_num: TInt32, iter_num: TInt32) -> TArray(TFloat, 1):
        """
        # todo: document
        Arguments:
            mode_num    (int32)         : the modulation mode.
            iter_num    (int32)         : the optimization iteration.
        Returns:
            # todo: document - list of abs angle
        """
        voltage_min_v, voltage_max_v = (0., 0.)
        if iter_num >= 100:
            # use error to set the scan range
            voltage_step_v = 0.
            if mode_num == 0:   voltage_step_v = self.dc_voltage_optima_0[iter_num, 2] * 10
            else:               voltage_step_v = self.dc_voltage_optima_1[iter_num, 2] * 10

            # set voltage scan bounds based on num_steps
            voltage_min_v = self.dc_voltage_optima_running[mode_num] - voltage_step_v * self.num_steps / 2.
            voltage_max_v = self.dc_voltage_optima_running[mode_num] + voltage_step_v * self.num_steps / 2.
        else:
            # set the voltage scan bounds based on voltage range
            voltage_min_v, voltage_max_v = self.dc_scan_range_volts_list[mode_num]

        # todo: adjust attenuation based on correlated amplitude

        # ensure voltages are within bounds
        voltage_min_v = max(self.dc_scan_range_volts_list[mode_num, 0], voltage_min_v)
        voltage_max_v = min(self.dc_scan_range_volts_list[mode_num, 1], voltage_max_v)
        print('\n\t\tVoltage range: [{:.2f}, {:.2f}]'.format(voltage_min_v, voltage_max_v))

        # create voltage scan array
        voltage_scan_arr = np.linspace(voltage_min_v, voltage_max_v, self.num_steps)
        np.random.shuffle(voltage_scan_arr)

        # reset intermediate result holders and related data structures
        self._res_holder_iter = 0
        self._res_holder_tmp = np.zeros((len(voltage_scan_arr) * self.repetitions_per_voltage, 3), dtype=float)

        return voltage_scan_arr

    @kernel(flags={"fast-math"})
    def _scan_voltage_mode(self, mode_freq_ftw: TInt32, mode_att_mu: TInt32, mode_dc_channel_num: TInt32,
                           voltage_scan_v_arr: TArray(TFloat, 1)):
        """
        Characterize micromotion along a mode.
        Arguments:
            mode_freq_ftw       (int32)             : the modulation frequency (as a 32-bit frequency tuning word).
            mode_att_my         (int32)             : the attenuation to use (in machine units).
            mode_dc_channel_num (int32)             : the voltage channel number to be scanned.
            voltage_scan_v_arr  (TArray(TFloat, 1)) : the list of voltages to scan.
        """
        self.core.break_realtime()

        # prepare modulation DDS
        self.dds_parametric.set_att_mu(mode_att_mu)
        self.dds_parametric.set_mu(mode_freq_ftw, asf=self.dds_parametric.ampl_modulation_asf,
                                   profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.core.break_realtime()

        # iterate over repetitions per voltage
        for rep_num in range(self.repetitions_per_voltage):

            # scan voltage configurations in the voltage vector
            for voltage_v in voltage_scan_v_arr:

                # set DC voltage
                self.voltage_set(mode_dc_channel_num, voltage_v)
                # synchronize hardware clock with timeline, then add delay for voltages to settle
                self.core.wait_until_mu(now_mu())
                delay_mu(self.time_dc_synchronize_delay_mu)

                # run parametric excitation and get timestamps
                pmt_timestamp_list = self.parametric_subsequence.run()

                # demodulate counts and store
                self.core.break_realtime()
                self._demodulate_counts(mode_freq_ftw, voltage_v, pmt_timestamp_list)
                self.core.reset()

            # rescue ion as needed
            self.rescue_subsequence.run(rep_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

    @rpc(flags={"async"})
    def _demodulate_counts(self, freq_mu: TInt32, voltage_v: TFloat, timestamp_mu_list: TArray(TInt64, 1)):
        """
        Demodulate the PMT timestamps with respect to the parametric modulation frequency.
        Arguments:
            freq_ftw            (int32)             : the modulation frequency (as a 32-bit frequency tuning word).
            voltage_v           (float)             : the current shim voltage (in volts).
            timestamp_mu_list   (TArray(TInt64, 1)) : the list of timestamps (in machine units) to demodulate.
        """
        # convert arguments to more convenient units
        freq_mhz = self.dds_parametric.ftw_to_frequency(freq_mu) / MHz
        timestamps_s = self.core.mu_to_seconds(np.array(timestamp_mu_list))
        
        # digitally demodulate counts and convert correlated signal to polar coordinates
        correlated_signal = np.mean(np.exp((2.j * np.pi * freq_mhz * 1e6) * timestamps_s))
        correlated_ampl =   np.abs(correlated_signal)
        correlated_phase =  np.angle(correlated_signal)
        # extract count rate in seconds
        count_rate_hz =     len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # store results in global parent dataset
        self.update_results(freq_mhz, voltage_v, correlated_ampl, correlated_phase, count_rate_hz)
        # store results relevant to optimization in an intermediate dataset
        self._res_holder_tmp[self._res_holder_iter] = np.array([voltage_v, correlated_ampl, correlated_phase])
        self._res_holder_iter += 1

    @rpc(flags={"async"})
    def _extract_optimum(self, res_arr: TArray(TFloat, 2)) -> TArray(TFloat, 1):
        """
        # todo: document
        Arguments:
            res_arr (int32): todo: document
        Returns:
            (float, float): the extracted voltage optimum and its error.
        """
        # group array by voltage and average demodulated signal
        res_arr_tmp = res_arr[np.argsort(res_arr[:, 0])]
        res_arr_averaged = []
        for key, group in groupby(res_arr_tmp, lambda row: row[0]):
            grouped_result = np.array(list(group))
            res_arr_averaged.append(np.mean(grouped_result, axis=0))

        # format results into a 2D array with complex type for complex linear fitting
        res_arr_averaged = np.array(res_arr_averaged)
        results_tmp = np.array([
            res_arr_averaged[:, 0],
            res_arr_averaged[:, 1] * np.exp(1.j * res_arr_averaged[:, 2])
        ], dtype='complex128').transpose()
        # sort results by voltage
        results_tmp = results_tmp[np.argsort(results_tmp[:, 0], axis=0)]

        # extract minimum mode voltage
        opt_voltage_v, opt_voltage_err = complexLinearFitMinimize(results_tmp)
        return opt_voltage_v, opt_voltage_err

    def analyze(self):
        pass

