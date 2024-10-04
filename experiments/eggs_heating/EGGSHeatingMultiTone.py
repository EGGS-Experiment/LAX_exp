import numpy as np
from artiq.experiment import *

from collections import OrderedDict
from itertools import product, chain, zip_longest

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)
import pandas as pd


class EGGSHeatingMultiTone(LAXExperiment, Experiment):
    """
    Experiment: EGGS Heating Multi Tone

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    Scans multiple ranges in an incredible way and allows user choice of amplitudes.
    """
    name = 'EGGS Heating'

    OSCILLATOR_CHANNELS = 5

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=40, ndecimals=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config", BooleanValue(default=True))
        self.setattr_argument("sub_repetitions", NumberValue(default=1, ndecimals=0, step=1, min=1, max=500))

        # # get subsequences
        self.initialize_subsequence = InitializeQubit(self)
        self.sidebandcool_subsequence = SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence = SidebandReadout(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # EGGS RF
        self.setattr_argument("scan_configuration", EnumerationValue(["Normal", "Randomize"], default="Normal"),
                              group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_mhz_list_0", Scannable(
            default=[
                RangeScan(81, 82, 300, randomize=False),
                ExplicitScan([3.]),
                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_mhz_list_1", Scannable(
            default=[
                RangeScan(83, 84, 300, randomize=False),
                ExplicitScan([3.]),
                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_mhz_list_2", Scannable(
            default=[
                RangeScan(85, 86, 300, randomize=False),
                ExplicitScan([3.]),
                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')

        self.setattr_argument("freq_eggs_heating_mhz_list_3", Scannable(
            default=[
                RangeScan(87, 88, 300, randomize=False),
                ExplicitScan([3.]),
                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')

        self.setattr_argument("freq_eggs_heating_mhz_list_4", Scannable(
            default=[
                RangeScan(90, 91, 300, randomize=False),
                ExplicitScan([3.]),
                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')

        self.setattr_argument("freq_eggs_heating_secular_khz_list", Scannable(
            default=[
                ExplicitScan([777.5]),
                CenterScan(777.5, 4, 0.1, randomize=True),
                ExplicitScan([767.2, 319.2, 1582, 3182]),
            ],
            global_min=0, global_max=10000, global_step=1,
            unit="kHz", scale=1, ndecimals=3
        ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_readout_us_list", Scannable(
            default=[
                ExplicitScan([98.8]),
                RangeScan(0, 1500, 100, randomize=True), ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, ndecimals=5
        ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("time_eggs_heating_ms",
                              NumberValue(default=1.0, ndecimals=5, step=1, min=0.000001, max=100000),
                              group='EGGS_Heating.waveform.time_phase')

        self.setattr_argument("phase_eggs_heating_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 9, randomize=True),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, ndecimals=3
        ), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("enable_amplitude_calibration", BooleanValue(default=True),
                              group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("att_eggs_heating_db", NumberValue(default=20., ndecimals=1, step=0.5, min=0, max=31.5),
                              group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_tone_pct_0",
                              NumberValue(default=12., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')

        self.setattr_argument("ampl_eggs_heating_tone_pct_1",
                              NumberValue(default=12., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')

        self.setattr_argument("ampl_eggs_heating_tone_pct_2",
                              NumberValue(default=12., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')

        self.setattr_argument("ampl_eggs_heating_tone_pct_3",
                              NumberValue(default=12., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')

        self.setattr_argument("ampl_eggs_heating_tone_pct_4",
                              NumberValue(default=12., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'),
                              group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",
                              NumberValue(default=100, ndecimals=1, step=100, min=10, max=100000),
                              group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",
                              NumberValue(default=500, ndecimals=0, step=100, min=100, max=2000),
                              group='EGGS_Heating.pulse_shaping')

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')
        #
        # # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')
        # tmp remove

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        self.num_configs = 1

        print(self.ampl_eggs_heating_tone_pct_1)

        if (self.ampl_eggs_heating_tone_pct_0 + self.ampl_eggs_heating_tone_pct_1 + self.ampl_eggs_heating_tone_pct_2 +
            self.ampl_eggs_heating_tone_pct_3 + self.ampl_eggs_heating_tone_pct_4) > 99.:
            raise Exception("Total Amplitude Should not Exceed 99%")

        max_freq_length = np.max((len(self.freq_eggs_heating_mhz_list_0), len(self.freq_eggs_heating_mhz_list_1),
                                  len(self.freq_eggs_heating_mhz_list_2),
                                  len(self.freq_eggs_heating_mhz_list_3), len(self.freq_eggs_heating_mhz_list_4)))

        '''SUBSEQUENCE PARAMETERS'''

        # convert build arguements into numpy arrays for use later
        freq_eggs_heating_mhz_list_0 = np.array(list(self.freq_eggs_heating_mhz_list_0))
        freq_eggs_heating_mhz_list_1 = np.array(list(self.freq_eggs_heating_mhz_list_1))
        freq_eggs_heating_mhz_list_2 = np.array(list(self.freq_eggs_heating_mhz_list_2))
        freq_eggs_heating_mhz_list_3 = np.array(list(self.freq_eggs_heating_mhz_list_3))
        freq_eggs_heating_mhz_list_4 = np.array(list(self.freq_eggs_heating_mhz_list_4))

        # get readout values
        self.freq_sideband_readout_ftw_list = np.array([100, 110])
        # self.freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.freq_sideband_readout_ftw_list = np.array(list(self.freq_sideband_readout_ftw_list))
        self.time_readout_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_readout_us_list])

        '''EGGS HEATING - TIMING'''
        self.time_eggs_heating_mu = self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser sample period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        if self.time_eggs_heating_mu % self.phaser_eggs.t_sample_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_sample_multiples = round(self.time_eggs_heating_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_eggs_heating_mu = np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        # add delay time after EGGS pulse to allow RF servo to re-lock
        self.time_rf_servo_holdoff_mu = self.get_parameter("time_rf_servo_holdoff_us", group="eggs",
                                                           conversion_function=us_to_mu)

        '''EGGS HEATING - PHASES'''
        # preallocate variables for phase
        self.phase_ch1_turns = np.float64(0)
        self.phase_phaser_turns_arr = np.zeros((2, 5), dtype=float)

        '''EGGS HEATING - CONFIG'''

        # convert build arguments to appropriate values and format as numpy arrays
        self.freq_eggs_secular_hz_list = np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz
        self.phase_eggs_heating_turns_list = np.array(list(self.phase_eggs_heating_turns_list))

        # create config data structure with amplitude values
        self.config_eggs_heating_list = np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                  max_freq_length *
                                                  len(self.freq_eggs_secular_hz_list) *
                                                  len(self.phase_eggs_heating_turns_list) *
                                                  len(self.time_readout_mu_list),
                                                  14), dtype=float)

        # pad frequency ranges that are shorter than the longest frequency range with zeros to ensure equal length
        self.freq_eggs_heating_mhz_list_0_padded = np.append(freq_eggs_heating_mhz_list_0,
                                                             np.zeros(
                                                                 max_freq_length - len(freq_eggs_heating_mhz_list_0)))
        self.freq_eggs_heating_mhz_list_1_padded = np.append(freq_eggs_heating_mhz_list_1,
                                                             np.zeros(
                                                                 max_freq_length - len(freq_eggs_heating_mhz_list_1)))
        self.freq_eggs_heating_mhz_list_2_padded = np.append(freq_eggs_heating_mhz_list_2,
                                                             np.zeros(
                                                                 max_freq_length - len(freq_eggs_heating_mhz_list_2)))
        self.freq_eggs_heating_mhz_list_3_padded = np.append(freq_eggs_heating_mhz_list_3,
                                                             np.zeros(
                                                                 max_freq_length - len(freq_eggs_heating_mhz_list_3)))
        self.freq_eggs_heating_mhz_list_4_padded = np.append(freq_eggs_heating_mhz_list_4,
                                                             np.zeros(
                                                                 max_freq_length - len(freq_eggs_heating_mhz_list_4)))

        # create an amplitude list is zero if we have gone outside one of the frequency ranges
        self.ampl_eggs_heating_tone_pct_0_list = (self.freq_eggs_heating_mhz_list_0_padded > 0) * self.ampl_eggs_heating_tone_pct_0
        self.ampl_eggs_heating_tone_pct_1_list = (self.freq_eggs_heating_mhz_list_1_padded > 0) * self.ampl_eggs_heating_tone_pct_1
        self.ampl_eggs_heating_tone_pct_2_list = (self.freq_eggs_heating_mhz_list_2_padded > 0) * self.ampl_eggs_heating_tone_pct_2
        self.ampl_eggs_heating_tone_pct_3_list = (self.freq_eggs_heating_mhz_list_3_padded > 0) * self.ampl_eggs_heating_tone_pct_3
        self.ampl_eggs_heating_tone_pct_4_list = (self.freq_eggs_heating_mhz_list_4_padded > 0) * self.ampl_eggs_heating_tone_pct_4

        # set padded values equal to the last value in the original list (avoid setting phaser freq to 0)
        self.freq_eggs_heating_mhz_list_0_padded[len(freq_eggs_heating_mhz_list_0):] = freq_eggs_heating_mhz_list_0[-1]
        self.freq_eggs_heating_mhz_list_1_padded[len(freq_eggs_heating_mhz_list_1):] = freq_eggs_heating_mhz_list_1[-1]
        self.freq_eggs_heating_mhz_list_2_padded[len(freq_eggs_heating_mhz_list_2):] = freq_eggs_heating_mhz_list_2[-1]
        self.freq_eggs_heating_mhz_list_3_padded[len(freq_eggs_heating_mhz_list_3):] = freq_eggs_heating_mhz_list_3[-1]
        self.freq_eggs_heating_mhz_list_4_padded[len(freq_eggs_heating_mhz_list_4):] = freq_eggs_heating_mhz_list_4[-1]

        # get the frequencies the oscillator of the phaser will be set to
        self.freq_eggs_heating_hz_list_0_padded = np.array(list(self.freq_eggs_heating_mhz_list_0_padded)) * MHz
        self.diff_freqs_hz_1 = (self.freq_eggs_heating_mhz_list_1_padded - self.freq_eggs_heating_mhz_list_0_padded) * MHz
        self.diff_freqs_hz_2 = (self.freq_eggs_heating_mhz_list_2_padded - self.freq_eggs_heating_mhz_list_0_padded) * MHz
        self.diff_freqs_hz_3 = (self.freq_eggs_heating_mhz_list_3_padded - self.freq_eggs_heating_mhz_list_0_padded) * MHz
        self.diff_freqs_hz_4 = (self.freq_eggs_heating_mhz_list_4_padded - self.freq_eggs_heating_mhz_list_0_padded) * MHz

        if (np.abs(np.array([self.diff_freqs_hz_1, self.diff_freqs_hz_2,
                             self.diff_freqs_hz_3, self.diff_freqs_hz_4])) > 10 * MHz).any():
            raise Exception("Ensure ALL Frequencies are within a 10 MHz window")

        # start populating eggs config list
        self.config_eggs_heating_list[:, [0, 10, 11, 12, 13]] = np.stack(
            np.meshgrid(self.freq_eggs_heating_hz_list_0_padded,
                        self.freq_eggs_secular_hz_list,
                        self.phase_eggs_heating_turns_list,
                        self.time_readout_mu_list,
                        self.freq_sideband_readout_ftw_list),
            -1).reshape(-1, 5)

        # ensure these frequency arrays are the same shape as the grid created above
        for i, arr in enumerate(np.array([self.diff_freqs_hz_1, self.diff_freqs_hz_2,
                                          self.diff_freqs_hz_3, self.diff_freqs_hz_4])):
            arr_reshaped = np.stack(np.meshgrid(arr, self.freq_eggs_secular_hz_list,
                                                self.phase_eggs_heating_turns_list,
                                                self.time_readout_mu_list,
                                                self.freq_sideband_readout_ftw_list), -1).reshape(-1, 5)
            self.config_eggs_heating_list[:, i + 1] = arr_reshaped[:, 0]

        # ensure these amplitude arrays are the same shape as the grid created above
        for i, arr in enumerate(np.array([self.ampl_eggs_heating_tone_pct_0_list, self.ampl_eggs_heating_tone_pct_1_list,
                          self.ampl_eggs_heating_tone_pct_1_list, self.ampl_eggs_heating_tone_pct_3_list,
                          self.ampl_eggs_heating_tone_pct_2_list])):
            arr_reshaped = np.stack(np.meshgrid(arr, self.freq_eggs_secular_hz_list,
                                                self.phase_eggs_heating_turns_list,
                                                self.time_readout_mu_list,
                                                self.freq_sideband_readout_ftw_list), -1).reshape(-1, 5)
            self.config_eggs_heating_list[:, i + 5] = arr_reshaped[:, 0]

        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config:                               np.random.shuffle(self.config_eggs_heating_list)

        '''EGGS HEATING - AMPLITUDE CALIBRATION'''
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points = self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve = Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # TMP REMOVE: MAKE SURE SIDEBAND AMPLITUDES ARE SCALED CORRECTLY FOLLOWING USER INPUT SPECS
        # todo: move to a phaser internal function
        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            pass

        '''EGGS HEATING - EGGS RF CONFIGURATION'''
        # configure pulse shaping
        # note: instead of having to deal with adjusting shape, etc., will just add the pulse shaping in addition to the actual pulse
        self._prepare_pulseshape()

        # configure for no phase-shift keying
        self.phaser_run = self.phaser_run_nopsk

        self.set_dataset("config", self.config_eggs_heating_list)
        # print(self.config_eggs_heating_list)
        # raise Exception("stop here")

    def _prepare_pulseshape(self):
        """
        Calculate waveform and timings for phaser pulse shaping.
        """
        ### PULSE SHAPING - TIMING ###
        # convert build variables to units of choice
        self.time_pulse_shape_rolloff_mu = self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)

        # todo: document
        self.t_max_phaser_update_rate_mu = 25 * self.phaser_eggs.t_sample_mu
        self.time_pulse_shape_sample_mu = self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
        # todo: add error handling and printouts if there's some sampling problem

        # ensure pulse shaping time is a multiple of the max sustained phaser update rate
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators) is (conservatively) about 1.5 MSPS (i.e. 25 sample periods))
        if self.time_pulse_shape_sample_mu % self.t_max_phaser_update_rate_mu:
            # round pulse shaping sample time up to the nearest multiple of phaser sample period
            t_sample_multiples = round((self.time_pulse_shape_sample_mu / self.t_max_phaser_update_rate_mu) + 0.5)
            self.time_pulse_shape_sample_mu = np.int64(t_sample_multiples * self.t_max_phaser_update_rate_mu)

        # delay time between successive updates to the pulse envelope, accounts for 2x t_sample_mu delay from having to set 3 oscillators
        self.time_pulse_shape_delay_mu = np.int64(self.time_pulse_shape_sample_mu - 2 * self.phaser_eggs.t_sample_mu)
        # note: calculation of the number of samples accounts for the delay from setting multiple oscillators
        self.num_pulse_shape_samples = np.int32(self.time_pulse_shape_rolloff_mu / (self.time_pulse_shape_sample_mu))

        # create holder object for pulse amplitudes
        self.ampl_pulse_shape_frac_list = np.tile(np.array([self.ampl_eggs_heating_tone_pct_0,
                                                            self.ampl_eggs_heating_tone_pct_1,
                                                            self.ampl_eggs_heating_tone_pct_2,
                                                            self.ampl_eggs_heating_tone_pct_3,
                                                            self.ampl_eggs_heating_tone_pct_4, ]) / 100.,
                                                  self.num_pulse_shape_samples).reshape(-1, 5)

        # calculate windowing values
        if self.type_pulse_shape == 'sine_squared':
            # calculate sine squared window
            self.ampl_window_frac_list = np.power(np.sin(
                (np.pi / (2. * self.num_pulse_shape_samples)) *
                np.linspace(1, self.num_pulse_shape_samples, self.num_pulse_shape_samples)),
                2)
        elif self.type_pulse_shape == 'error_function':
            raise Exception('Error: error function window not implemented')
        else:
            raise Exception('Error: idk, some window problem')

        # todo: add other spicy windows
        # ensure window array has correct dimensions required
        self.ampl_window_frac_list = np.array([self.ampl_window_frac_list]).transpose()

        # apply window to pulse shape
        self.ampl_pulse_shape_frac_list *= self.ampl_window_frac_list
        self.ampl_pulse_shape_reverse_frac_list = self.ampl_pulse_shape_frac_list[::-1]

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions * len(self.config_eggs_heating_list),
                6)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        ### PHASER INITIALIZATION ###
        self.phaser_setup()
        self.core.break_realtime()

        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)

        # reset debug triggers
        self.ttl8.off()
        self.ttl9.off()
        # tmp remove
        self.ttl10.off()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_eggs_pulseshape_rise = self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        _handle_eggs_pulseshape_fall = self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # tmp remove
        # set phaser attenuators
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)
        self.core.break_realtime()
        # tmp remove

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # implement sub-repetitions here to avoid initial overhead
            _subrep_iter = 0
            _config_iter = 0

            # sweep experiment configurations
            while _config_iter < self.num_configs:

                for freq_readout_ftw in self.freq_sideband_readout_ftw_list:
                    '''CONFIGURE'''
                    config_vals = self.config_eggs_heating_list[_config_iter]
                    # extract values from config list
                    sideband_freq_hz = config_vals[10]
                    phase_turns = config_vals[11]
                    time_readout_mu = np.int64(config_vals[12])
                    freq_readout_ftw = np.int32(config_vals[13])
                    ampl_tone_0 = config_vals[5]
                    ampl_tone_1 = config_vals[6]
                    ampl_tone_2 = config_vals[7]
                    ampl_tone_3 = config_vals[8]
                    ampl_tone_4 = config_vals[9]

                    center_freq = config_vals[0]
                    diff_freq_1 = config_vals[1]
                    diff_freq_2 = config_vals[2]
                    diff_freq_3 = config_vals[3]
                    diff_freq_4 = config_vals[4]

                    self.core.break_realtime()
                    # configure EGGS tones and set readout frequency
                    self.phaser_configure(center_freq, diff_freq_1, diff_freq_2, diff_freq_3, diff_freq_4, phase_turns)
                    self.core.break_realtime()
                    self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                      profile=0)
                    self.core.break_realtime()

                    '''STATE PREPARATION'''
                    # initialize ion in S-1/2 state
                    self.initialize_subsequence.run_dma()
                    # sideband cool
                    self.sidebandcool_subsequence.run_dma()

                    '''EGGS HEATING'''
                    # EGGS - START/SETUP
                    # activate integrator hold
                    self.ttl10.on()
                    # # set phaser attenuators
                    # at_mu(self.phaser_eggs.get_next_frame_mu())
                    # self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
                    # delay_mu(self.phaser_eggs.t_sample_mu)
                    # self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

                    # reset DUC phase to start DUC deterministically
                    self.phaser_eggs.reset_duc_phase()
                    self.core_dma.playback_handle(_handle_eggs_pulseshape_rise)

                    # EGGS - RUN
                    self.phaser_run(ampl_tone_0, ampl_tone_1, ampl_tone_2, ampl_tone_3, ampl_tone_4)

                    # EGGS - STOP
                    self.core_dma.playback_handle(_handle_eggs_pulseshape_fall)
                    self.phaser_eggs.phaser_stop()
                    # deactivate integrator hold
                    self.ttl10.off()
                    # add delay time after EGGS pulse to allow RF servo to re-lock
                    delay_mu(self.time_rf_servo_holdoff_mu)

                    '''READOUT'''
                    self.sidebandreadout_subsequence.run_time(time_readout_mu)
                    self.readout_subsequence.run_dma()
                    counts = self.readout_subsequence.fetch_count()

                    # update dataset
                    with parallel:
                        self.update_results(
                            freq_readout_ftw,
                            counts,
                            center_freq,
                            sideband_freq_hz,
                            phase_turns,
                            time_readout_mu
                        )
                        self.core.break_realtime()

                    # resuscitate ion
                    self.rescue_subsequence.resuscitate()

                    # death detection
                    self.rescue_subsequence.detect_death(counts)
                    self.core.break_realtime()

                # check termination more frequently in case reps are low
                if (_loop_iter % 50) == 0:
                    self.check_termination()
                    self.core.break_realtime()
                _loop_iter += 1

                # handle sub-repetition logic
                if _config_iter % 2 == 1:
                    _subrep_iter += 1
                    if _subrep_iter < self.sub_repetitions:
                        _config_iter -= 1
                    else:
                        _subrep_iter = 0
                        _config_iter += 1
                else:
                    _config_iter += 1

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        '''CLEANUP'''
        self.core.break_realtime()
        self.phaser_eggs.reset_oscillators()
        # tmp remove
        self.ttl10.off()
        # tmp remove

    '''
    HELPER FUNCTIONS - PHASER
    '''

    @kernel(flags={"fast-math"})
    def phaser_setup(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        phase_turns = self.phase_eggs_heating_turns
        center_freq = self.config_eggs_heating_list[0, 0]
        diff_freq_1 = self.config_eggs_heating_list[0, 1]
        diff_freq_2 = self.config_eggs_heating_list[0, 2]
        diff_freq_3 = self.config_eggs_heating_list[0, 3]
        diff_freq_4 = self.config_eggs_heating_list[0, 4]
        self.core.break_realtime()
        # configure EGGS tones and set readout frequency
        self.phaser_configure(center_freq, diff_freq_1, diff_freq_2, diff_freq_3, diff_freq_4, phase_turns)

        # record phaser rising pulse shape DMA sequence
        self.core.break_realtime()
        if self.enable_pulse_shaping:
            with self.core_dma.record('_PHASER_PULSESHAPE_RISE'):
                # set amplitude values at given time
                for ampl_val_list in self.ampl_pulse_shape_frac_list:
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2],
                                                 ampl_val_list[3], ampl_val_list[4])
                    delay_mu(self.time_pulse_shape_delay_mu)
        else:
            with self.core_dma.record('_PHASER_PULSESHAPE_RISE'):
                pass

        # record phaser falling pulse shape DMA sequence
        self.core.break_realtime()
        if self.enable_pulse_shaping:
            with self.core_dma.record('_PHASER_PULSESHAPE_FALL'):
                # set amplitude values at given time
                for ampl_val_list in self.ampl_pulse_shape_reverse_frac_list:
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2],
                                                 ampl_val_list[3], ampl_val_list[4])
                    delay_mu(self.time_pulse_shape_delay_mu)
        else:
            with self.core_dma.record('_PHASER_PULSESHAPE_FALL'):
                pass

    @kernel(flags={"fast-math"})
    def phaser_configure(self, center_freq_hz: TFloat, diff_freq_hz_1: TFloat, diff_freq_hz_2: TFloat,
                         diff_freq_hz_3: TFloat, diff_freq_hz_4: TFloat, phase_turns: TFloat) -> TNone:
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            center_freq_hz         (float)     : the carrier frequency (in Hz).
            diff_freqs_hz      (float)     : the sideband frequency (in Hz).
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        self.phase_ch1_turns = (self.phaser_eggs.phase_inherent_ch1_turns +
                                (center_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0
        self.phase_phaser_turns_arr[0, 0] = phase_turns
        self.phase_phaser_turns_arr[1, 0] = phase_turns
        # oscillator 1
        self.phase_phaser_turns_arr[0, 1] = phase_turns
        self.phase_phaser_turns_arr[1, 1] = phase_turns
        # oscillator 2
        self.phase_phaser_turns_arr[0, 2] = phase_turns
        self.phase_phaser_turns_arr[1, 2] = phase_turns
        # oscillator 3
        self.phase_phaser_turns_arr[0, 3] = phase_turns
        self.phase_phaser_turns_arr[1, 3] = phase_turns
        # oscillator 4
        self.phase_phaser_turns_arr[0, 4] = phase_turns
        self.phase_phaser_turns_arr[1, 4] = phase_turns
        self.core.break_realtime()

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(center_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(center_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(self.phase_ch1_turns)
        # todo: do I need to add another get_next_frame_mu?
        self.phaser_eggs.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 - duc freq
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(diff_freq_hz_1)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(diff_freq_hz_1)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(diff_freq_hz_2)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(diff_freq_hz_2)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 3
        with parallel:
            self.phaser_eggs.channel[0].oscillator[3].set_frequency(diff_freq_hz_3)
            self.phaser_eggs.channel[1].oscillator[3].set_frequency(diff_freq_hz_3)
            delay_mu(self.phaser_eggs.t_sample_mu)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 4
        with parallel:
            self.phaser_eggs.channel[0].oscillator[4].set_frequency(diff_freq_hz_4)
            self.phaser_eggs.channel[1].oscillator[4].set_frequency(diff_freq_hz_4)
            delay_mu(self.phaser_eggs.t_sample_mu)
            delay_mu(self.phaser_eggs.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_pulseshape_point(self, ampl_tone_0_pct: TFloat, ampl_tone_1_pct: TFloat, ampl_tone_2_pct: TFloat,
                                ampl_tone_3_pct: TFloat, ampl_tone_4_pct: TFloat) -> TNone:
        """
        Convenience function for pulse-shaping the EGGS waveform.
        Sets the same oscillator amplitudes for both channels.
        The phase used is the existing value of self.phaser_phaser_turns_arr.
        Arguments:
            ampl_tone_0_pct (float) : amplitude of oscillator 0 (as a decimal fraction).
            ampl_tone_1_pct (float) : amplitude of oscillator 1 (as a decimal fraction).
            ampl_tone_2_pct (float) : amplitude of oscillator 2 (as a decimal fraction).
            ampl_tone_3_pct (float) : amplitude of oscillator 3 (as a decimal fraction).
            ampl_tone_4_pct (float) : amplitude of oscillator 4 (as a decimal fraction).
        """
        # set oscillator 0
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_tone_0_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 0],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_tone_0_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 0],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_tone_1_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 1],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_tone_1_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 1],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_tone_2_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_tone_2_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)
        # set oscillator 3
        with parallel:
            self.phaser_eggs.channel[0].oscillator[3].set_amplitude_phase(amplitude=ampl_tone_3_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[3].set_amplitude_phase(amplitude=ampl_tone_3_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)
        # set oscillator 4
        with parallel:
            self.phaser_eggs.channel[0].oscillator[4].set_amplitude_phase(amplitude=ampl_tone_4_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[4].set_amplitude_phase(amplitude=ampl_tone_4_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)

    '''
    HELPER FUNCTIONS - PSK
    '''

    @kernel(flags={"fast-math"})
    def phaser_run_nopsk(self, ampl_tone_0_pct: TFloat, ampl_tone_1_pct: TFloat, ampl_tone_2_pct: TFloat,
                         ampl_tone_3_pct: TFloat, ampl_tone_4_pct: TFloat) -> TNone:
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same oscillator amplitudes for both channels.
        Arguments:
            ampl_tone_0_pct (float) : amplitude of oscillator 0 (as a decimal fraction).
            ampl_tone_1_pct (float) : amplitude of oscillator 1 (as a decimal fraction).
            ampl_tone_2_pct (float) : amplitude of oscillator 2 (as a decimal fraction).
            ampl_tone_3_pct (float) : amplitude of oscillator 3 (as a decimal fraction).
            ampl_tone_4_pct (float) : amplitude of oscillator 4 (as a decimal fraction).
        """
        # set oscillator 0
        with parallel:
            self.ttl8.on()
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_tone_0_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 0],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_tone_0_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 0],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_tone_1_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 1],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_tone_2_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 1],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_tone_2_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_tone_2_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)
        # set oscillator 3
        with parallel:
            self.phaser_eggs.channel[0].oscillator[3].set_amplitude_phase(amplitude=ampl_tone_3_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[3].set_amplitude_phase(amplitude=ampl_tone_3_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)
        # set oscillator 4
        with parallel:
            self.phaser_eggs.channel[0].oscillator[4].set_amplitude_phase(amplitude=ampl_tone_4_pct,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[4].set_amplitude_phase(amplitude=ampl_tone_4_pct,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)

        # main eggs pulse
        delay_mu(self.time_eggs_heating_mu)
        self.ttl8.off()

    '''
    ANALYSIS
    '''

    def analyze(self):
        # should be backward-compatible with experiment which had no sub-repetitions
        try:
            sub_reps = self.sub_repetitions
        except Exception as e:
            sub_reps = 1

        # handle errors from data processing
        try:
            # print results
            print("\tResults - EGGS Heating:")

            # convert dataset to array
            dataset = np.array(self.results)

            ## determine scan type
            # carrier sweep
            if len(np.unique(dataset[:, 2])) > 1:
                sorting_col_num = 2
            # secular frequency
            elif len(np.unique(dataset[:, 3])) > 1:
                sorting_col_num = 3
            else:
                sorting_col_num = 0

            ## convert results to sideband ratio
            ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
                                                                                           1, 0,
                                                                                           self.repetitions, sub_reps)
            phonons = convert_ratios_to_coherent_phonons(ratios)

            ## process secular frequency sweep
            if sorting_col_num == 3:
                fit_params_secular, fit_err_secular, _ = fitSincGeneric(scanning_freq_MHz, phonons)
                phonon_n = fit_params_secular[0]
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_secular', fit_params_secular)
                self.set_dataset('fit_err_secular', fit_err_secular)

                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_secular, fit_err_secular]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False,
                                 archive=False)

                # print results to log
                print(
                    "\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_secular[1] * 1e3, fit_err_secular[1] * 1e3))


            ## process sideband readout sweep
            else:
                rsb_freqs_MHz, bsb_freqs_MHz, _ = extract_sidebands_freqs(scanning_freq_MHz)
                fit_params_rsb, fit_err_rsb, fit_rsb = fitSincGeneric(rsb_freqs_MHz, ave_rsb)
                fit_params_bsb, fit_err_bsb, fit_bsb = fitSincGeneric(bsb_freqs_MHz, ave_bsb)
                phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_rsb', fit_params_rsb)
                self.set_dataset('fit_params_bsb', fit_params_bsb)
                self.set_dataset('fit_err_rsb', fit_err_rsb)
                self.set_dataset('fit_err_bsb', fit_err_bsb)

                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False,
                                 archive=False)

                # print results to log
                print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2., float(fit_err_rsb[1]) / 2.))
                print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2., float(fit_err_bsb[1]) / 2.))

        except Exception as e:
            print("Warning: unable to process data.")
            print(repr(e))