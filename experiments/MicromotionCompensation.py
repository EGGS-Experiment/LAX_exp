import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.system.subsequences import ParametricExcite
import LAX_exp.experiments.ParametricSweep as ParametricSweep

import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config


class MicromotionCompensation(ParametricSweep.ParametricSweep, Experiment):
    """
    Experiment: Micromotion Compensation

    # todo: redocument
    Modulate the trap RF close to a secular frequency while sweeping shim voltgaes
    to measure micromotion.
    """
    name = 'Micromotion Compensation'


    def build_experiment(self):
        # get DC channel configuration dictionary
        self.dc_config_channeldict =                                dc_config.channeldict

        # core arguments
        self.setattr_argument("iterations",                         NumberValue(default=2, ndecimals=0, step=1, min=1, max=100))
        self.setattr_argument("adaptive",                           BooleanValue(default=False))


        # modulation - general
        self.setattr_argument("repetitions_per_voltage",            NumberValue(default=5, ndecimals=0, step=1, min=1, max=100), group='modulation_general')
        self.setattr_argument("att_mod_db",                         NumberValue(default=18, ndecimals=1, step=0.5, min=0, max=31.5), group='modulation_general')

        # modulation - mode #1
        self.setattr_argument("freq_mod0_khz",                      NumberValue(default=1315, ndecimals=3, step=10, min=1, max=10000), group='modulation_1')
        self.setattr_argument("dc_channel_mod0",                    EnumerationValue(list(self.dc_config_channeldict.keys()), default='V Shim'), group='modulation_1')
        self.setattr_argument("dc_voltages_mod0_v_list",            Scannable(
                                                                            default=[
                                                                                CenterScan(45., 20., 2., randomize=True),
                                                                                ExplicitScan([40.])
                                                                            ],
                                                                            global_min=0, global_max=200, global_step=1,
                                                                            unit="V", scale=1, ndecimals=1
                                                                        ), group='modulation_1')
        # self.setattr_argument("dc_voltages_mod0_v_range",           PYONValue([10, 80]), group='modulation_1')
        # self.setattr_argument("dc_voltages_mod0_v_step",        NumberValue(default=1055, ndecimals=3, step=10, min=1, max=10000), group='modulation_1')
        # self.setattr_argument("dc_voltages_mod0_num_points",        NumberValue(default=1055, ndecimals=3, step=10, min=1, max=10000), group='modulation_1')

        # modulation - mode #2
        self.setattr_argument("freq_mod1_khz",                      NumberValue(default=1100, ndecimals=3, step=10, min=1, max=10000), group='modulation_2')
        self.setattr_argument("dc_channel_mod1",                    EnumerationValue(list(self.dc_config_channeldict.keys()), default='H Shim'), group='modulation_2')
        self.setattr_argument("dc_voltages_mod1_v_list",            Scannable(
                                                                        default=[
                                                                            CenterScan(45., 20., 2., randomize=True),
                                                                            ExplicitScan([40.])
                                                                        ],
                                                                        global_min=0, global_max=200, global_step=1,
                                                                        unit="V", scale=1, ndecimals=1
                                                                    ), group='modulation_2')

        # cooling
        self.setattr_argument("ampl_cooling_pct",                   NumberValue(default=50, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",                   NumberValue(default=110, ndecimals=6, step=1, min=1, max=500), group='cooling')

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_modulation')

        # get relevant subsequences
        self.parametric_subsequence =                               ParametricExcite(self)

    def prepare_experiment(self):
        # tmp remove
        # todo: get min, max, step, num_points
        # self.dc_voltages_mod0_v_list

        # convert cooling parameters to machine units
        self.freq_cooling_ftw =                                     self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.ampl_cooling_asf =                                     self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.time_cooling_holdoff_mu =                              self.core.seconds_to_mu(3 * ms)

        # convert general modulation parameters to machine units
        self.att_modulation_mu =                                    att_to_mu(self.att_mod_db * dB)

        # convert modulation (mode #1) parameters to machine units
        self.freq_mod0_ftw =                                        self.dds_modulation.frequency_to_ftw(self.freq_mod0_khz * kHz)
        self.dc_channel_mod0_num =                                  self.dc_config_channeldict[self.dc_channel_mod0]['num']
        self.dc_voltages_mod0_v_list =                              np.array(list(self.dc_voltages_mod0_v_list))

        # convert modulation (mode #2) parameters to machine units
        self.freq_mod1_ftw =                                        self.dds_modulation.frequency_to_ftw(self.freq_mod1_khz * kHz)
        self.dc_channel_mod1_num =                                  self.dc_config_channeldict[self.dc_channel_mod1]['num']
        self.dc_voltages_mod1_v_list =                              np.array(list(self.dc_voltages_mod1_v_list))

        # create array to store current voltage values
        # note: use median of voltage scan range for starting voltages
        self.dc_voltage_optimal_v_list =                            np.array([np.median(self.dc_voltages_mod0_v_list),
                                                                              np.median(self.dc_voltages_mod1_v_list)])
        # create vector for storing voltage basis direction vector
        # note: mostly used for adaptive optimization
        self.dc_voltage_basis_vector =                              np.array([[1., 0.], [0., 1.]], dtype=float)

        # tmp remove
        # create holder dataset
        max_size_tmp =                                              max(len(self.dc_voltages_mod0_v_list), len(self.dc_voltages_mod1_v_list))
        self.tmp_dataset_holder =                                   np.zeros((self.repetitions_per_voltage * max_size_tmp, 4), dtype=float)
        self.tmp_dataset_iter =                                     0
        self._tmp_result_iter =                                     0

        # connect to labrad and get DC server
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                   self.cxn.dc_server

    @property
    def results_shape(self):
        return (self.iterations * self.repetitions_per_voltage * (len(self.dc_voltages_mod0_v_list) + len(self.dc_voltages_mod1_v_list)),
                6)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # run given number of iterations
        for _iter_num in range(self.iterations):
            self.core.break_realtime()

            # generate voltage vector array (mode 0) and optimize
            voltage_0_list_tmp = self._generate_voltage_vector_array(0)
            self.core.break_realtime()
            self._sweep_voltage_mode(self.freq_mod0_ftw, 0, voltage_0_list_tmp)
            self.core.break_realtime()
            # todo: put process here?
            # todo: get resultant voltage/save idk

            # generate voltage vector array (mode 1) and optimize
            voltage_1_list_tmp = self._generate_voltage_vector_array(1)
            self.core.break_realtime()
            self._sweep_voltage_mode(self.freq_mod1_ftw, 1, voltage_1_list_tmp)
            # todo: put process here?
            # todo: get resultant voltage/save idk

            # todo: guess/adapt voltage vector basis
            if self.adaptive:
                pass


    # ANALYZE
    def analyze(self):
        # todo: fit data for all voltage values
        # todo: extract uncertainties
        # todo: if multiple voltages, return optimal voltage
        pass


    # HELPER FUNCTIONS
    @rpc
    def _generate_voltage_vector_array(self, mode_num: TInt32) -> TArray(TFloat, 2):
        """
        # todo: document

        Arguments:
            mode_num            (int)         : the modulation mode number tmp remove idk
        Returns:
            # todo: document - list of abs angle
        """
        # todo: get current values?
        # todo: fancier return function
        # self.dc_voltage_basis_vector[mode_num]
        array_return_tmp = np.array([])
        if mode_num == 0:
            array_return_tmp = np.zeros((len(self.dc_voltages_mod0_v_list), 2), dtype=float)
            array_return_tmp[:, 0] = self.dc_voltages_mod0_v_list
            array_return_tmp[:, 1] = self.dc_voltage_optimal_v_list[1]
        elif mode_num == 1:
            array_return_tmp = np.zeros((len(self.dc_voltages_mod1_v_list), 2), dtype=float)
            array_return_tmp[:, 1] = self.dc_voltages_mod1_v_list
            array_return_tmp[:, 0] = self.dc_voltage_optimal_v_list[0]

        return array_return_tmp

    @kernel(flags={"fast-math"})
    def _sweep_voltage_mode(self, mode_freq_ftw: TInt32, mode_num: TInt32,
                            voltage_vec_v_arr: TArray(TFloat, 2)) -> TArray(TFloat, 1):
        """
        todo: document

        Arguments:
            mode_freq_ftw            (int32)         : the modulation frequency (as a 32-bit frequency tuning word).
            idk tmp remove      (float)         : the current voltage for modulation channel 2 (in volts).
            voltage_vec_v_arr      (float)         : the current voltage for modulation channel 1 (in volts).
        Returns:
            # todo: document - list of abs angle
        """
        # set modulation frequency for the given mode
        self.dds_modulation.set_mu(mode_freq_ftw, asf=self.dds_modulation.ampl_modulation_asf)
        self.core.break_realtime()

        # iterate over repetitions per voltage
        for rep_num in range(self.repetitions_per_voltage):

            # sweep voltage configurations in the voltage vector
            for voltage_vec_v in voltage_vec_v_arr:

                # extract voltage values from voltage vector
                voltage_mod0_v = voltage_vec_v[0]
                voltage_mod1_v = voltage_vec_v[1]
                self.core.break_realtime()

                # set mode voltages
                self.voltage_set(self.dc_channel_mod0_num, voltage_mod0_v)
                self.voltage_set(self.dc_channel_mod1_num, voltage_mod1_v)
                self.core.break_realtime()

                # add holdoff period for recooling the ion
                delay_mu(self.time_cooling_holdoff_mu)
                # todo: see if processing can be made parallel with holdoff time (end)

                # run parametric excitation and get timestamps
                pmt_timestamp_list = self.parametric_subsequence.run()

                # demodulate counts and store them in the global results holder, as well as a local data structure
                with parallel:
                    self.core.reset()
                    self._process_results(mode_freq_ftw, voltage_mod0_v, voltage_mod1_v, pmt_timestamp_list)

        # do complex linear fit to extract optimal voltage
        opt_voltages_v_list = self._extract_optimum(mode_num)
        return opt_voltages_v_list

    @rpc
    def _process_results(self, freq_mu: TInt32,
                         voltage_mod0_v: TFloat, voltage_mod1_v: TFloat,
                         timestamp_mu_list: TArray(TInt64, 1)):
        """
        Convert modulation frequency and timestamps from machine units and demodulate.

        Arguments:
            freq_ftw            (int32)         : the modulation frequency (as a 32-bit frequency tuning word).
            voltage_mod0_v      (float)         : the current voltage for modulation channel 1 (in volts).
            voltage_mod1_v      (float)         : the current voltage for modulation channel 2 (in volts).
            timestamp_mu_list   (array(int64))  : the list of timestamps (in machine units) to demodulate.
        Returns:
            # todo: document - list of abs angle
        """
        # convert frequency from ftw to mhz
        freq_mhz = self.dds_modulation.ftw_to_frequency(freq_mu) / MHz

        # convert timestamps and digitally demodulate counts
        timestamps_s = self.core.mu_to_seconds(np.array(timestamp_mu_list))
        correlated_signal = np.mean(np.exp((2.j * np.pi * freq_mhz * 1e6) * timestamps_s))

        # convert demodulated signal to polar coordinates (i.e. abs and angle)
        correlated_ampl = np.abs(correlated_signal)
        correlated_phase = np.angle(correlated_signal)
        # get count rate in seconds
        count_rate_hz = len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # update dataset and return result
        self.update_results(freq_mhz, voltage_mod0_v, voltage_mod1_v, correlated_ampl, correlated_phase, count_rate_hz)

        # also store values in a local dataset
        self.tmp_dataset_holder[self.tmp_dataset_iter] = np.array([voltage_mod0_v, voltage_mod1_v, correlated_ampl, correlated_phase])
        self.tmp_dataset_iter += 1

    @rpc
    def _extract_optimum(self, mode_num: TInt32) -> TArray(TFloat, 1):
        """
        # todo: document
        # todo: do we need mode as arg?

        Arguments:
            mode_num            (int32)         : the modulation frequency (as a 32-bit frequency tuning word).
        Returns:
            # todo: document - list of abs angle
        """
        # todo: document
        results_tmp = np.array(self.tmp_dataset_holder[: self.tmp_dataset_iter], dtype='float')

        # todo: groupBy and average for given mode
        _reduce_func = lambda data: np.array([np.mean(data[:, 1]), np.mean(data[:, 2])])
        results_tmp = groupBy(results_tmp, column_num=mode_num, reduce_func=_reduce_func)
        results_tmp = np.concatenate((np.array([list(results_tmp.keys())]).transpose(),
                                      np.array(list(results_tmp.values()))),
                                     axis=-1)

        # format results into a 2D array with complex type for complex linear fitting
        results_tmp = np.array([
            results_tmp[:, 0],
            results_tmp[:, 1] * np.exp(1.j * results_tmp[:, 2])
        ], dtype='complex128').transpose()

        # sanitize data holder objects
        self.tmp_dataset_iter = 0
        self.tmp_dataset_holder *= 0

        # tmp remove
        # print(np.abs(results_tmp))
        self.set_dataset('iter{:d}_mode{}'.format(self._tmp_result_iter//2, mode_num), results_tmp)
        self._tmp_result_iter += 1
        # tmp remove

        # extract minimum mode voltage
        opt_voltage_v = complexFitMinimize(results_tmp)
        opt_voltage_list_v_tmp = self.dc_voltage_optimal_v_list

        # ensure values are within acceptable range, otherwise return original values
        if mode_num == 0:
            if (opt_voltage_v > np.max(self.dc_voltages_mod0_v_list)) | (opt_voltage_v < np.min(self.dc_voltages_mod0_v_list)):
                print('ERROR: MODE 0 OUTSIDE MAX VOLTAGE ({:.2f})'.format(opt_voltage_v))
                opt_voltage_v = self.dc_voltage_optimal_v_list[mode_num]

            opt_voltage_list_v_tmp[mode_num] = opt_voltage_v
            # todo: get the voltage for the other mode
            # opt_voltage_list_v_tmp[1] = _interpolateData_tmp(opt_voltage_v)
        elif mode_num == 1:
            if (opt_voltage_v > np.max(self.dc_voltages_mod1_v_list)) | (opt_voltage_v < np.min(self.dc_voltages_mod1_v_list)):
                print('ERROR: MODE 1 OUTSIDE MAX VOLTAGE ({:.2f})'.format(opt_voltage_v))
                opt_voltage_v = self.dc_voltage_optimal_v_list[mode_num]

            opt_voltage_list_v_tmp[mode_num] = opt_voltage_v
            # todo: get the voltage for the other mode
            # opt_voltage_list_v_tmp[0] = _interpolateData_tmp(opt_voltage_v)

        # set new optimal values and print result
        self.dc_voltage_optimal_v_list = opt_voltage_list_v_tmp
        print('\t\tOptimal Voltages: {:.2f}, {:.2f}'.format(opt_voltage_list_v_tmp[0], opt_voltage_list_v_tmp[1]))
        return opt_voltage_list_v_tmp

    # @rpc
    # def _extract_optimum(self, mode_num: TInt32, dataset: TArray(TFloat, 2)) -> TArray(TFloat, 1):
    #     """
    #     # todo: document
    #     # todo: do we need mode as arg?
    #
    #     Arguments:
    #         mode_num            (int32)         : the modulation frequency (as a 32-bit frequency tuning word).
    #         dataaset   (array(int64))  : the list of timestamps (in machine units) to demodulate.
    #     Returns:
    #         # todo: document - list of abs angle
    #     """
    #     # todo: groupBy and average for given mode
    #     results_tmp = np.array(dataset, dtype='float')
    #     results_tmp = groupBy(ds_tmp, column_num=mode_num, reduce_func=np.mean)
    #     results_tmp = np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()
    #
    #     # format results into a 2D array with complex type for complex linear fitting
    #     ds_tmp = np.array(dataset, dtype='complex128')
    #     # todo: more elegant way of formatting?
    #     if mode_num == 0:
    #         results_tmp = np.array([
    #             results_tmp[:, 0],
    #             results_tmp[:, 2] * np.exp(1.j * results_tmp[:, 3])
    #         ], dtype='complex').transpose()
    #     elif mode_num== 1:
    #         results_tmp = np.array([
    #             results_tmp[:, 1],
    #             results_tmp[:, 2] * np.exp(1.j * results_tmp[:, 3])
    #         ], dtype='complex').transpose()
    #
    #     # todo: interpolate the other mode_num and get what the voltage should be
    #     m0, b0 = (0., 0.)
    #     points = ds_tmp[[0, -1]]
    #     if mode_num == 0:
    #         m0 = (points[1, 1] - points[1, 0]) / (points[0, 1] - points[0, 0])
    #         b0 = points[0, 1] - m0 * points[0, 0]
    #     elif mode_num == 1:
    #         m0 = (points[0, 1] - points[0, 0]) / (points[0, 1] - points[0, 0])
    #         b0 = points[0, 1] - m0 * points[0, 0]
    #
    #     # create interpolation function
    #     def _interpolateData_tmp(x):
    #         """
    #         todo: document
    #         """
    #         return m0 * x + b0
    #
    #
    #     # extract minimum mode voltage
    #     opt_voltage_v = complexFitMinimize(ds_tmp)
    #     opt_voltage_list_v_tmp = np.array([0., 0.], dtype=float)
    #
    #     # todo: ensure values are within acceptable range, otherwise return original values
    #     if mode_num == 0:
    #         if (opt_voltage_v > np.max(self.dc_voltages_mod0_v_list)) | (opt_voltage_v < np.min(self.dc_voltages_mod0_v_list)):
    #             print('ERROR: MODE 0 OUTSIDE MAX VOLTAGE ({:.2f})'.format(opt_voltage_v))
    #             opt_voltage_v = self.dc_voltage_optimal_v_list[mode_num]
    #
    #         opt_voltage_list_v_tmp[mode_num] = opt_voltage_v
    #         # todo: get the voltage for the other mode
    #         opt_voltage_list_v_tmp[1] = _interpolateData_tmp(opt_voltage_v)
    #     elif mode_num == 1:
    #         if (opt_voltage_v > np.max(self.dc_voltages_mod1_v_list)) | (opt_voltage_v < np.min(self.dc_voltages_mod1_v_list)):
    #             print('ERROR: MODE 1 OUTSIDE MAX VOLTAGE ({:.2f})'.format(opt_voltage_v))
    #             opt_voltage_v = self.dc_voltage_optimal_v_list[mode_num]
    #
    #         opt_voltage_list_v_tmp[mode_num] = opt_voltage_v
    #         # todo: get the voltage for the other mode
    #         opt_voltage_list_v_tmp[0] = _interpolateData_tmp(opt_voltage_v)
    #
    #     # todo: print result
    #     print('\t\tOptimal Voltages: {:.2f}, {:.2f}'.format(opt_voltage_list_v_tmp[0], opt_voltage_list_v_tmp[1]))
    #     return opt_voltage_list_v_tmp
