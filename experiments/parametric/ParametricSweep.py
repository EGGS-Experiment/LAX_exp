import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config

import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import ParametricExcite, RescueIon


class InsufficientCounts(Exception):
    """
    Raised when the PMT is not getting sufficient counts for some reason.
    Prevents fixed PMT count-type experiments (e.g. ParametricSweep) from
    blocking due to insufficient PMT counts.
    """
    pass


class ParametricSweep(LAXExperiment, Experiment):
    """
    Experiment: Parametric Sweep

    Modulate the trap RF close to a secular frequency while sweeping shim voltgaes
    to measure micromotion.
    """
    name = 'Parametric Sweep'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=15, ndecimals=1, step=0.5, min=0, max=31.5), group='modulation')
        self.setattr_argument("mod_freq_khz_list",                  Scannable(
                                                                        default= [
                                                                            ExplicitScan([1252.4]),
                                                                            CenterScan(1090.8, 10, 0.25, randomize=True),
                                                                        ],
                                                                        global_min=1, global_max=200000, global_step=1,
                                                                        unit="kHz", scale=1, ndecimals=4
                                                                    ), group='modulation')

        # shimming voltages
        self.dc_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_channel",             EnumerationValue(list(self.dc_channeldict.keys()), default='H Shim'), group='voltage')
        self.setattr_argument("dc_voltages_v_list",     Scannable(
                                                                        default=[
                                                                            CenterScan(52.8, 4., 0.1, randomize=True),
                                                                            ExplicitScan([38.5]),
                                                                        ],
                                                                        global_min=0, global_max=400, global_step=1,
                                                                        unit="V", scale=1, ndecimals=1
                                                                    ), group='voltage')

        # cooling
        self.setattr_argument("ampl_cooling_pct",                   NumberValue(default=25, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",                   NumberValue(default=105, ndecimals=6, step=1, min=1, max=500), group='cooling')

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_parametric')
        self.setattr_device('pmt')

        # get relevant subsequences
        self.parametric_subsequence =                               ParametricExcite(self)
        self.rescue_subsequence =                                   RescueIon(self)


    def prepare_experiment(self):
        # get voltage parameters
        self.dc_channel_num =                               self.dc_channeldict[self.dc_channel]['num']
        self.dc_voltages_v_list =                           np.array(list(self.dc_voltages_v_list))
        self.time_dc_synchronize_delay_mu =                 self.core.seconds_to_mu(488 * ms)

        # convert cooling parameters to machine units
        self.ampl_cooling_asf =                             self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.freq_cooling_ftw =                             self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.time_cooling_holdoff_mu =                      self.core.seconds_to_mu(3 * ms)

        # modulation control and synchronization
        self.att_modulation_mu =                            att_to_mu(self.mod_att_db * dB)
        self.freq_modulation_list_mu =                      np.array([
                                                                self.dds_parametric.frequency_to_ftw(freq_mhz * kHz)
                                                                for freq_mhz in self.mod_freq_khz_list
                                                            ])

        # connect to labrad
        self.cxn =                                          labrad.connect(environ['LABRADHOST'],
                                                                           port=7682, tls_mode='off',
                                                                           username='', password='lab')
        self.dc =                                           self.cxn.dc_server

        # set up variables for ensuring PMT counts are above some threshold
        self.fluorescence_calibration_time_mu =             np.int64(30000000)  # 30ms
        self.fluorescence_calibration_threshold_counts =    400


    @property
    def results_shape(self):
        return (self.repetitions * len(self.dc_voltages_v_list) * len(self.mod_freq_khz_list),
                5)


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel: TInt32, voltage_v: TFloat) -> TNone:
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        voltage_set_v = self.dc.voltage_fast(channel, voltage_v)
        # print('\tvoltage set: {}'.format(voltage_set_v))

    @rpc(flags={"async"})
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)
        self.dc.serial_write('remote.w 1\r\n')
        self.dc.serial_read('\n')


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
            with sequential:
                self.dds_parametric.set_att_mu(self.att_modulation_mu)
                self.dds_parametric.set_phase_absolute()

            # set up labrad devices via RPC
            self.prepareDevicesLabrad()
        self.core.break_realtime()

        # do check to verify that mirror is flipped to mirror
        # tmp remove: fix
        # self._check_PMT_counting()
        # tmp remove: fix

    @kernel(flags={"fast-math"})
    def _check_PMT_counting(self):
        """
        Check that everything is set up to ensure the ion is fluorescing correctly.
        Runs a given number of PMT readout sequences and compares the average fluorescence
        against some given number.
        """
        # count fluorescence
        self.pmt.count(self.fluorescence_calibration_time_mu)

        # ensure fluorescence exceeds threshold
        counts_calibration = self.pmt.fetch_count()
        if counts_calibration < self.fluorescence_calibration_threshold_counts:
            raise InsufficientCounts("Error: PMT not receiving sufficient counts.")


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # run given number of repetitions
        for trial_num in range(self.repetitions):

            # sweep voltage
            for voltage_v in self.dc_voltages_v_list:

                # set DC voltage
                self.voltage_set(self.dc_channel_num, voltage_v)

                # synchronize hardware clock with timeline, then add delay for voltages to settle
                # note: delay has added advantage of recooling the ion
                self.core.wait_until_mu(now_mu())
                delay_mu(self.time_dc_synchronize_delay_mu)

                # sweep modulation frequencies
                for freq_mu in self.freq_modulation_list_mu:
                    self.core.break_realtime()

                    # set modulation frequency
                    self.dds_parametric.set_mu(freq_mu, asf=self.dds_parametric.ampl_modulation_asf,
                                               profile=0, phase_mode=PHASE_MODE_CONTINUOUS)

                    # run parametric excitation
                    pmt_timestamp_list = self.parametric_subsequence.run()

                    # process results (stores them in our results dataset for us)
                    with parallel:
                        self._process_results(freq_mu,
                                              voltage_v,
                                              pmt_timestamp_list)
                        self.core.reset()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()


    @rpc(flags={"async"})
    def _process_results(self, freq_mu: TInt32, voltage_v: TFloat, timestamp_mu_list: TArray(TInt64, 1)):
        """
        Convert modulation frequency and timestamps from machine units and demodulate.

        Arguments:
            freq_ftw            (int32)         : the modulation frequency (as a 32-bit frequency tuning word).
            voltage_v           (float)         : the current shim voltage (in volts).
            timestamp_mu_list   (list(int64))   : the list of timestamps (in machine units) to demodulate.
        """
        # convert frequency to mhz
        freq_mhz = self.dds_parametric.ftw_to_frequency(freq_mu) / MHz

        # convert timestamps and digitally demodulate counts
        timestamps_s = self.core.mu_to_seconds(np.array(timestamp_mu_list))
        correlated_signal = np.mean(np.exp((2.j * np.pi * freq_mhz * 1e6) * timestamps_s))
        # convert demodulated signal to polar coordinates (i.e. abs and angle)
        correlated_ampl = np.abs(correlated_signal)
        correlated_phase = np.angle(correlated_signal)

        # extract count rate in seconds
        count_rate_hz = len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # update dataset
        self.update_results(freq_mhz, voltage_v, correlated_ampl, correlated_phase, count_rate_hz)


    # ANALYSIS
    def analyze_experiment(self):
        """
        Fit resultant spectra with a sinc profile to extract n,
        then fit a line to extract heating rate
        """
        # todo: group first by freq and voltage, and mean
        # todo: then fit for each freq (account for disjunction if both eggs and rf somehow)
        # todo: then fit for each freq
        # note: no need to convert units since we already do this in _process_results
        # also, we leave out counts (5th column) from fitting
        results_tmp =           np.array(self.results)[:, :4]
        # group dataset first by shim voltage
        results_tmp =           groupBy(results_tmp, column_num=1)

        # group dataset by modulation frequency and average results,
        # then convert resultant dict into 2D dataset
        _reduce_func =          lambda ds: np.mean(ds, axis=0)
        _dict_to_array =        lambda _dict: np.concatenate((np.array([list(_dict.keys())]).transpose(),
                                                              np.array(list(_dict.values()))),
                                                             axis=-1)
        results_tmp =           {key_voltage: _dict_to_array(groupBy(val_dataset, column_num=0, reduce_func=_reduce_func))
                                 for key_voltage, val_dataset in results_tmp.items()}

        # extract optimal voltage for each modulation frequency
        if len(results_tmp) > 1:
            from itertools import groupby
            def groupByList(dataset, column_num=0, reduce_func=lambda x: x):
                # sort array by given column number (necessary for itertools.groupby)
                dataset = np.array(dataset)
                dataset = dataset[np.argsort(dataset[:, column_num])]

                # group same values into dict
                return [
                    [key, reduce_func(np.array([val for val in group]))]
                    for key, group in groupby(dataset, lambda arr: arr[column_num])
                ]

            # group results by frequency, then sort by voltage
            __results_tmp = np.array([sublist[1] for sublist in groupByList(np.array(self.results), column_num=1)])[::-1]
            __results_tmp = np.array([sublist[np.argsort(sublist[:, 0])] for sublist in __results_tmp])

            def extractVoltageOptimum(col_num):
                # create array of [voltage, a*exp(i \theta)]
                _data = __results_tmp[:, col_num, [1, 2, 3]]
                _data = np.array([
                    _data[:, 0],
                    _data[:, 1] * np.exp(1.j * _data[:, 2])
                ], dtype='complex128').transpose()
                return complexLinearFitMinimize(_data)

            # process voltage optima for all frequencies
            voltage_optima = [
                [__results_tmp[0, col_num, 0] * 1000.,  *extractVoltageOptimum(col_num)]
                for col_num in range(len(__results_tmp[0, :, 0]))
            ]
            # np.array(list(voltage_optima.items()))

            # save and print results
            self.set_dataset('voltage_optima_khz_v', voltage_optima)
            for res_arr in voltage_optima:
                print('\tFreq. (kHz): {:.2f}\tOpt. Voltage (V): {:.2f} +/- {:.2f} V'.format(res_arr[0], res_arr[1], res_arr[2]))

        # todo: check that frequencies are continuous
        if len(results_tmp) == 1:
            # fit amplitude for all voltages
            results_amplitude_fit = {key_voltage: fitDampedDrivenOscillatorAmplitude(val_dataset[:, [0, 1]])
                                     for key_voltage, val_dataset in results_tmp.items()}
            # todo: use amplitude fit values to support phase fitting
            # results_phase_fit =     {key_voltage: fitDampedDrivenOscillatorPhase(val_dataset[:, [0, 2]])
            #                          for key_voltage, val_dataset in results_tmp.items()}

            # process fit results
            amplitude_fit_params_saved =   np.array([[key_voltage, *fit_params[0]]
                                                     for key_voltage, fit_params in results_amplitude_fit.items()])
            amplitude_fit_err_saved =       np.array([[key_voltage, *fit_params[1]]
                                                      for key_voltage, fit_params in results_amplitude_fit.items()])

            # save results to hdf5 as a dataset
            self.set_dataset('amplitude_fit_params_saved',  amplitude_fit_params_saved)
            self.set_dataset('amplitude_fit_err_saved',     amplitude_fit_err_saved)

            # save results to dataset manager for dynamic experiments
            res_dj = [amplitude_fit_params_saved, amplitude_fit_err_saved]
            self.set_dataset('temp.parametricsweep.results', res_dj, broadcast=True, persist=False, archive=False)
            self.set_dataset('temp.parametricsweep.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

            # print out fit results
            print("\tResults - Parametric Sweep:")
            # todo: make it always print out the mean of results (or not? maybe?)
            if len(results_tmp.keys()) == 1:
                print("\t\tMode Frequency (kHz):\t{:.2f} +/- {:.2f}".format(amplitude_fit_params_saved[0, 2]*1.e3, amplitude_fit_err_saved[0, 2]*1.e3))
                print("\t\tLinewidth (kHz):\t{:.2f} +/- {:.2f}".format(abs(amplitude_fit_params_saved[0, 3])*1.e3, amplitude_fit_err_saved[0, 3]*1.e3))
            else:
                # todo
                pass
