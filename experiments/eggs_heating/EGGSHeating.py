import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)


class EGGSHeating(LAXExperiment, Experiment):
    """
    Experiment: EGGS Heating

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'EGGS Heating'
    kernel_invariants = {
        # config-related
        'freq_sideband_readout_ftw_list', 'time_readout_mu_list', 'freq_eggs_carrier_hz_list', 'freq_eggs_secular_hz_list',
        'phase_eggs_heating_rsb_turns_list',
        'time_eggs_heating_mu', 'time_rf_servo_holdoff_mu',
        'config_eggs_heating_list', 'num_configs',
        # pulse shaping
        'time_pulse_shape_rolloff_mu', 't_max_phaser_update_rate_mu', 'time_pulse_shape_sample_mu', 'time_pulse_shape_delay_mu',
        'num_pulse_shape_samples', 'ampl_pulse_shape_frac_list', 'ampl_window_frac_list', 'ampl_pulse_shape_frac_list', 'ampl_pulse_shape_reverse_frac_list',
        # PSK
        'config_dynamical_decoupling_psk_list', 'time_psk_delay_mu', 'phaser_run'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",        NumberValue(default=40, ndecimals=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",   BooleanValue(default=True))
        self.setattr_argument("sub_repetitions",    NumberValue(default=1, ndecimals=0, step=1, min=1, max=500))

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandcool_subsequence =     SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # EGGS RF
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    RangeScan(84.65, 84.95, 100, randomize=False),
                                                                                    ExplicitScan([3.]),
                                                                                    ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                                                                                    CenterScan(83.20175, 0.05, 0.0005, randomize=True),
                                                                                ],
                                                                                global_min=0.005, global_max=4800, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([1113.81]),
                                                                                    CenterScan(777.5, 4, 0.1, randomize=True),
                                                                                    ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_readout_us_list",                       Scannable(
                                                                                default=[
                                                                                    ExplicitScan([127.1]),
                                                                                    RangeScan(0, 1500, 100, randomize=True),
                                                                                ],
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, ndecimals=5
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=1.0, ndecimals=5, step=1, min=0.000001, max=100000), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 9, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",               NumberValue(default=0., ndecimals=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=14., ndecimals=1, step=0.5, min=0, max=31.5), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",                  NumberValue(default=38., ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",                  NumberValue(default=43., ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        # EGGS RF - waveform - amplitude - dynamical decoupling - configuration
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=True), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=0.4, ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",                       BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",                           EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",                NumberValue(default=100, ndecimals=1, step=100, min=10, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",                NumberValue(default=500, ndecimals=0, step=100, min=100, max=2000), group='EGGS_Heating.pulse_shaping')

        # EGGS RF - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_dd_phase_shift_keying",               BooleanValue(default=True), group='EGGS_Heating.waveform.psk')
        self.setattr_argument("num_dynamical_decoupling_phase_shifts",      NumberValue(default=3, ndecimals=0, step=10, min=1, max=100), group='EGGS_Heating.waveform.psk')

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')
        # tmp remove

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        # todo: ensure phaser amplitudes sum to less than 100%

        '''SUBSEQUENCE PARAMETERS'''
        # get readout values
        self.freq_sideband_readout_ftw_list =   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.time_readout_mu_list =             np.array([self.core.seconds_to_mu(time_us * us)
                                                          for time_us in self.time_readout_us_list])

        '''EGGS HEATING - TIMING'''
        self.time_eggs_heating_mu = self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser sample period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        # todo: move to an internal phaser function
        if self.time_eggs_heating_mu % self.phaser_eggs.t_sample_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_sample_multiples =        round(self.time_eggs_heating_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_eggs_heating_mu = np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        # add delay time after EGGS pulse to allow RF servo to re-lock
        self.time_rf_servo_holdoff_mu = self.get_parameter("time_rf_servo_holdoff_us", group="eggs",
                                                           conversion_function=us_to_mu)

        '''EGGS HEATING - PHASES'''
        # preallocate variables for phase
        self.phase_ch1_turns =          np.float(0)
        self.phase_phaser_turns_arr =   np.zeros((2, 3), dtype=float)


        '''EGGS HEATING - CONFIG'''
        # convert build arguments to appropriate values and format as numpy arrays
        self.freq_eggs_carrier_hz_list =            np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =            np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz
        self.phase_eggs_heating_rsb_turns_list =    np.array(list(self.phase_eggs_heating_rsb_turns_list))

        # create config data structure with amplitude values
        self.config_eggs_heating_list = np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                  len(self.freq_eggs_carrier_hz_list) *
                                                  len(self.freq_eggs_secular_hz_list) *
                                                  len(self.phase_eggs_heating_rsb_turns_list) *
                                                  len(self.time_readout_mu_list),
                                                  8), dtype=float)
        # note: sideband readout frequencies are at the end of the
        # meshgrid to support adjacent_sidebands configuration option
        self.config_eggs_heating_list[:, [1, 2, -2, -1, 0]] =   np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list,
                                                                                     self.freq_eggs_secular_hz_list,
                                                                                     self.phase_eggs_heating_rsb_turns_list,
                                                                                     self.time_readout_mu_list,
                                                                                     self.freq_sideband_readout_ftw_list),
                                                                         -1).reshape(-1, 5)
        self.config_eggs_heating_list[:, [3, 4, 5]] =           np.array([self.ampl_eggs_heating_rsb_pct,
                                                                          self.ampl_eggs_heating_bsb_pct,
                                                                          self.ampl_eggs_dynamical_decoupling_pct]) / 100.
        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config: np.random.shuffle(self.config_eggs_heating_list)

        # precalculate length of configuration list here to reduce run-time overhead
        self.num_configs = len(self.config_eggs_heating_list)


        '''EGGS HEATING - AMPLITUDE CALIBRATION'''
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                     self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                      Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # TMP REMOVE: MAKE SURE SIDEBAND AMPLITUDES ARE SCALED CORRECTLY FOLLOWING USER INPUT SPECS
        # todo: move to a phaser internal function
        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            for i, (_, carrier_freq_hz, secular_freq_hz, _, _, _, _, _) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz =                    (np.array([-secular_freq_hz, secular_freq_hz]) + carrier_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac =                        ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac
                # TMP FIX: MAKE SURE SCALED POWER FOLLOWS SPECIFICATIONS OF RSB AND BSB PCT
                # scaled_power_pct =                                          (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                #                                                              ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                scaled_power_pct =                              (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                                                                 ((self.ampl_eggs_heating_rsb_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                # update configs and convert amplitude to frac
                self.config_eggs_heating_list[i, [3, 4, 5]] =   np.array([scaled_power_pct[0],
                                                                          scaled_power_pct[1],
                                                                          self.ampl_eggs_dynamical_decoupling_pct]) / 100.


        '''EGGS HEATING - EGGS RF CONFIGURATION'''
        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling: self.config_eggs_heating_list[:, 5] = 0.

        # configure pulse shaping
        # note: instead of having to deal with adjusting shape, etc., will just add the pulse shaping in addition to the actual pulse
        self._prepare_pulseshape()

        # configure phase-shift keying for dynamical decoupling
        self._prepare_psk()

    def _prepare_pulseshape(self):
        """
        Calculate waveform and timings for phaser pulse shaping.
        """
        ### PULSE SHAPING - TIMING ###
        # convert build variables to units of choice
        self.time_pulse_shape_rolloff_mu =  self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)

        # todo: document
        self.t_max_phaser_update_rate_mu =  25 * self.phaser_eggs.t_sample_mu
        self.time_pulse_shape_sample_mu =   self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
        # todo: add error handling and printouts if there's some sampling problem

        # ensure pulse shaping time is a multiple of the max sustained phaser update rate
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators) is (conservatively) about 1.5 MSPS (i.e. 25 sample periods))
        if self.time_pulse_shape_sample_mu % self.t_max_phaser_update_rate_mu:
            # round pulse shaping sample time up to the nearest multiple of phaser sample period
            t_sample_multiples =                    round((self.time_pulse_shape_sample_mu / self.t_max_phaser_update_rate_mu) + 0.5)
            self.time_pulse_shape_sample_mu =       np.int64(t_sample_multiples * self.t_max_phaser_update_rate_mu)

        # delay time between successive updates to the pulse envelope, accounts for 2x t_sample_mu delay from having to set 3 oscillators
        self.time_pulse_shape_delay_mu =            np.int64(self.time_pulse_shape_sample_mu - 2 * self.phaser_eggs.t_sample_mu)
        # note: calculation of the number of samples accounts for the delay from setting multiple oscillators
        self.num_pulse_shape_samples =              np.int32(self.time_pulse_shape_rolloff_mu / (self.time_pulse_shape_sample_mu))

        ### PULSE SHAPING - AMPLITUDE WINDOW ###
        # hotfix - 2024/06/21 - implement dd enable for pulse shaping
        ampl_dd_pct =                               self.ampl_eggs_dynamical_decoupling_pct
        if not self.enable_dynamical_decoupling:    ampl_dd_pct = 0.

        # create holder object for pulse amplitudes
        self.ampl_pulse_shape_frac_list = np.tile(np.array([self.ampl_eggs_heating_rsb_pct,
                                                            self.ampl_eggs_heating_bsb_pct,
                                                            ampl_dd_pct]) / 100.,
                                                  self.num_pulse_shape_samples).reshape(-1, 3)

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
        self.ampl_pulse_shape_frac_list *=          self.ampl_window_frac_list
        self.ampl_pulse_shape_reverse_frac_list =   self.ampl_pulse_shape_frac_list[::-1]

        # create data structures to hold pulse shaping DMA sequences
        self.phaser_dma_handle_pulseshape_rise = (0, np.int64(0), np.int32(0))
        self.phaser_dma_handle_pulseshape_fall = (0, np.int64(0), np.int32(0))

    def _prepare_psk(self):
        """
        Calculate and prepare timings for PSK.
        """
        # create config holder for dynamical decoupling PSK; holds time_mu and phase in turns
        self.config_dynamical_decoupling_psk_list =         np.zeros((self.num_dynamical_decoupling_phase_shifts + 1, 2), dtype=np.int64)
        self.config_dynamical_decoupling_psk_list[1::2] =   1

        # divide total eggs heating time into PSK segments
        self.time_psk_delay_mu = np.int64(round(self.time_eggs_heating_mu / (self.num_dynamical_decoupling_phase_shifts + 1)))
        # ensure PSK interval time is a multiple of the phaser sample period
        if self.time_psk_delay_mu % self.phaser_eggs.t_sample_mu:
            # round dynamical decoupling PSK interval to the nearest multiple of phaser sample period
            t_sample_multiples =        round(self.time_psk_delay_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_psk_delay_mu =    np.int64(t_sample_multiples * self.phaser_eggs.t_sample_mu)

        # update dynamical decoupling config list with verified PSK time
        self.config_dynamical_decoupling_psk_list[:, 0] = self.time_psk_delay_mu
        # ensure that psk rate doesn't exceed the shaping time (t_max_phaser_update_rate_mu; about 25 * t_sample_mu)
        if self.enable_dd_phase_shift_keying and (self.time_psk_delay_mu < self.t_max_phaser_update_rate_mu):
            raise Exception("Error: num_dynamical_decoupling_phase_shifts too high; PSK update rate exceeds max sustained event rate.")

        # set appropriate phaser run method for dynamical decoupling PSK
        if self.enable_dd_phase_shift_keying:   self.phaser_run = self.phaser_run_psk
        else:                                   self.phaser_run = self.phaser_run_nopsk

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions * len(self.config_eggs_heating_list),
                6)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
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
    def run_main(self) -> TNone:
        self.core.break_realtime()

        # get custom sequence handles
        self.phaser_dma_handle_pulseshape_rise = self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        self.phaser_dma_handle_pulseshape_fall = self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # set phaser attenuators
        # note: this is done here instead of during sequence
        # since attenuator setting glitches cause heating if there is no
        # high-pass to filter them
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)
        self.core.break_realtime()


        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # implement sub-repetitions here to avoid initial overhead
            _subrep_iter = 0
            _config_iter = 0

            # sweep experiment configurations
            while _config_iter < self.num_configs:

                '''CONFIGURE'''
                config_vals = self.config_eggs_heating_list[_config_iter]
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                carrier_freq_hz =   config_vals[1]
                sideband_freq_hz =  config_vals[2]
                ampl_rsb_frac =     config_vals[3]
                ampl_bsb_frac =     config_vals[4]
                ampl_dd_frac =      config_vals[5]
                phase_rsb_turns =   config_vals[6]
                time_readout_mu =   np.int64(config_vals[7])
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz, phase_rsb_turns)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
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
                self.core_dma.playback_handle(self.phaser_dma_handle_pulseshape_rise)

                # EGGS - RUN
                self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

                # EGGS - STOP
                self.core_dma.playback_handle(self.phaser_dma_handle_pulseshape_fall)
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
                self.update_results(
                    freq_readout_ftw,
                    counts,
                    carrier_freq_hz,
                    sideband_freq_hz,
                    phase_rsb_turns,
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

            # support graceful termination at the repetition level
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
        # workaround: get starting phase values for pulse shaping
        carrier_freq_hz =   self.config_eggs_heating_list[0, 1]
        sideband_freq_hz =  self.config_eggs_heating_list[0, 2]
        phase_rsb_turns =   self.config_eggs_heating_list[0, 6]
        self.core.break_realtime()

        # configure EGGS tones and set readout frequency; also necessary to ensure phase delays are correctly set
        self.phaser_configure(carrier_freq_hz, sideband_freq_hz, phase_rsb_turns)

        # record phaser rising pulse shape DMA sequence
        self.core.break_realtime()
        if self.enable_pulse_shaping:
            with self.core_dma.record('_PHASER_PULSESHAPE_RISE'):
                # set amplitude values at given time
                for ampl_val_list in self.ampl_pulse_shape_frac_list:
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2])
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
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2])
                    delay_mu(self.time_pulse_shape_delay_mu)
        else:
            with self.core_dma.record('_PHASER_PULSESHAPE_FALL'):
                pass

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, phase_rsb_turns: TFloat) -> TNone:
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the carrier frequency (in Hz).
            sideband_freq_hz        (float)     : the sideband frequency (in Hz).
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        self.phase_ch1_turns =          (self.phaser_eggs.phase_inherent_ch1_turns +
                                         (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_phaser_turns_arr[0, 0] = phase_rsb_turns
        self.phase_phaser_turns_arr[1, 0] = phase_rsb_turns
        # oscillator 1 (BSB)
        self.phase_phaser_turns_arr[0, 1] = (sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        self.phase_phaser_turns_arr[1, 1] = (sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        # oscillator 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_phaser_turns_arr[0, 2] = 0.
        self.phase_phaser_turns_arr[1, 2] = 0.5
        self.core.break_realtime()

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
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
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_pulseshape_point(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat) -> TNone:
        """
        Convenience function for pulse-shaping the EGGS waveform.
        Sets the same oscillator amplitudes for both channels.
        The phase used is the existing value of self.phaser_phaser_turns_arr.
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[0, 0], clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[1, 0], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[0, 1], clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[1, 1], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[0, 2], clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[1, 2], clr=0)


    '''
    HELPER FUNCTIONS - PSK
    '''
    @kernel(flags={"fast-math"})
    def phaser_run_nopsk(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat) -> TNone:
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same oscillator amplitudes for both channels.
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        # set oscillator 0 (RSB)
        with parallel:
            self.ttl8.on()
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[0, 0], clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[1, 0], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[0, 1], clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[1, 1], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[0, 2], clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[1, 2], clr=0)

        # main eggs pulse
        delay_mu(self.time_eggs_heating_mu)
        self.ttl8.off()

    @kernel(flags={"fast-math"})
    def phaser_run_psk(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat) -> TNone:
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same RSB, BSB, and dynamical decoupling amplitudes for both channels.
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[0, 0], clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[1, 0], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[0, 1], clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[1, 1], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[0, 2], clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[1, 2], clr=0)

        # first PSK delay period
        delay_mu(self.config_dynamical_decoupling_psk_list[0][0])

        # conduct PSK on carrier
        for dd_config_vals in self.config_dynamical_decoupling_psk_list[1:]:
            # set oscillator 2 (carrier) with phase shift
            with parallel:
                self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                              phase=self.phase_phaser_turns_arr[0, 2] + (dd_config_vals[1] * 0.5),
                                                                              clr=0)
                self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                              phase=self.phase_phaser_turns_arr[1, 2] + (dd_config_vals[1] * 0.5),
                                                                              clr=0)
                delay_mu(dd_config_vals[0])


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
            if len(np.unique(dataset[:, 2])) > 1:       sorting_col_num = 2
            # secular frequency
            elif len(np.unique(dataset[:, 3])) > 1:     sorting_col_num = 3
            else:                                       sorting_col_num = 0

            ## convert results to sideband ratio
            ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num, 1, 0,
                                                                                           self.repetitions, sub_reps)
            phonons = convert_ratios_to_coherent_phonons(ratios)

            ## process secular frequency sweep
            if sorting_col_num == 3:
                fit_params_secular, fit_err_secular, _ = fitSincGeneric(scanning_freq_MHz, phonons)
                phonon_n = fit_params_secular[0]
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_secular',  fit_params_secular)
                self.set_dataset('fit_err_secular',     fit_err_secular)
                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_secular, fit_err_secular]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
                # print results to log
                print("\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_secular[1] * 1e3, fit_err_secular[1] * 1e3))

            ## process sideband readout sweep
            else:
                rsb_freqs_MHz, bsb_freqs_MHz, _ =       extract_sidebands_freqs(scanning_freq_MHz)
                fit_params_rsb, fit_err_rsb, fit_rsb =  fitSincGeneric(rsb_freqs_MHz, ave_rsb)
                fit_params_bsb, fit_err_bsb, fit_bsb =  fitSincGeneric(bsb_freqs_MHz, ave_bsb)
                phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_rsb',  fit_params_rsb)
                self.set_dataset('fit_params_bsb',  fit_params_bsb)
                self.set_dataset('fit_err_rsb',     fit_err_rsb)
                self.set_dataset('fit_err_bsb',     fit_err_bsb)
                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
                # print results to log
                print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2., float(fit_err_rsb[1]) / 2.))
                print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2., float(fit_err_bsb[1]) / 2.))

        except Exception as e:
            print("Warning: unable to process data.")
            print(repr(e))

