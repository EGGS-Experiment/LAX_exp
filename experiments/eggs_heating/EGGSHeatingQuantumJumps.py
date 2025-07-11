import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout,
                                         SidebandCoolContinuousRAM)
from sipyco import pyon


class EGGSHeatingQuantumJumps(LAXExperiment, Experiment):
    """
    Experiment: EGGS Heating Quantum Jumps

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'EGGS Heating Quantum Jumps'
    kernel_invariants = {
        # config-related
        'freq_sideband_readout_ftw_list', 'time_readout_mu_list', 'freq_eggs_carrier_hz_list', 'freq_eggs_secular_hz_list',
        'phase_eggs_heating_rsb_turns_list', 'time_eggs_heating_mu',
        'config_eggs_heating_list', 'num_configs',

        # pulse shaping
        'time_pulse_shape_rolloff_mu', 't_max_phaser_update_rate_mu', 'time_pulse_shape_sample_mu', 'time_pulse_shape_delay_mu',
        'num_pulse_shape_samples', 'ampl_pulse_shape_frac_list', 'ampl_window_frac_list', 'ampl_pulse_shape_frac_list', 'ampl_pulse_shape_reverse_frac_list',

        # PSK
        'config_dynamical_decoupling_psk_list', 'time_psk_delay_mu', 'phaser_run',

        # Quantum Jumps
        'freq_sideband_readout_mean_ftw', '_quantum_jump_monitor_rsb', '_quantum_jump_monitor_bsb',
        'config_quantum_jumps', '_results_quantum_jumps_idx', '_dataset_quantum_jumps_idx',
        'max_num_quantum_jumps', 'peak_sb_threshold_min_frac', 'peak_sb_threshold_max_frac',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",            NumberValue(default=1, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",       BooleanValue(default=True))
        self.setattr_argument("sub_repetitions",        NumberValue(default=40, precision=0, step=1, min=1, max=500))

        # quantum jump arguments
        self.setattr_argument("num_quantum_jumps",              NumberValue(default=500, precision=0, step=1, min=1, max=100000), group='EGGS_Heating.quantum_jumps')
        self.setattr_argument("count_threshold",                NumberValue(default=85, precision=0, step=10, min=0, max=250), group='EGGS_Heating.quantum_jumps')
        self.setattr_argument("quantum_jump_threshold_phonon",  NumberValue(default=0.75, precision=2, step=0.1, min=0., max=1.2), group='EGGS_Heating.quantum_jumps')

        # get subsequences
        self.initialize_subsequence =           InitializeQubit(self)
        # self.sidebandcool_subsequence =         SidebandCoolContinuous(self)

        self.profile_729_readout = 0
        self.profile_729_SBC = 1

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )
        self.sidebandreadout_subsequence = SidebandReadout(self, profile_dds=self.profile_729_readout)
        self.readout_subsequence =              Readout(self)
        self.rescue_subsequence =               RescueIon(self)

        # EGGS RF
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",     Scannable(
                                                                            default=[
                                                                                RangeScan(84.35, 84.55, 100, randomize=False),
                                                                                ExplicitScan([83.2028, 83.2028, 83.2028, 83.2028, 83.2097]),
                                                                                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
                                                                            ],
                                                                            global_min=0.005, global_max=4800, global_step=1,
                                                                            unit="MHz", scale=1, precision=6
                                                                        ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",     Scannable(
                                                                            default=[
                                                                                ExplicitScan([1117.16]),
                                                                                CenterScan(777.5, 4, 0.1, randomize=True),
                                                                                ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                            ],
                                                                            global_min=0, global_max=10000, global_step=1,
                                                                            unit="kHz", scale=1, precision=3
                                                                        ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_readout_us_list",                   Scannable(
                                                                            default=[
                                                                                ExplicitScan([123.3]),
                                                                                RangeScan(0, 1500, 100, randomize=True),
                                                                            ],
                                                                            global_min=1, global_max=100000, global_step=1,
                                                                            unit="us", scale=1, precision=5
                                                                        ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("time_eggs_heating_ms",                   NumberValue(default=1.0, precision=5, step=1, min=0.000001, max=100000), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns_list",      Scannable(
                                                                            default=[
                                                                                ExplicitScan([0.]),
                                                                                RangeScan(0, 1.0, 9, randomize=True),
                                                                            ],
                                                                            global_min=0.0, global_max=1.0, global_step=1,
                                                                            unit="turns", scale=1, precision=3
                                                                        ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",           NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("enable_amplitude_calibration",           BooleanValue(default=False), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("att_eggs_heating_db",                    NumberValue(default=14., precision=1, step=0.5, min=0, max=31.5), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",              NumberValue(default=31.5, precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",              NumberValue(default=43., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        # EGGS RF - waveform - amplitude - dynamical decoupling - configuration
        self.setattr_argument("enable_dynamical_decoupling",            BooleanValue(default=True), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",     NumberValue(default=0.4, precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=10, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=500, precision=0, step=100, min=100, max=2000), group='EGGS_Heating.pulse_shaping')

        # EGGS RF - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_dd_phase_shift_keying",           BooleanValue(default=True), group='EGGS_Heating.waveform.psk')
        self.setattr_argument("num_dynamical_decoupling_phase_shifts",  NumberValue(default=5, precision=0, step=10, min=1, max=100), group='EGGS_Heating.waveform.psk')

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
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
            t_sample_multiples =            round(self.time_eggs_heating_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_eggs_heating_mu =     np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        '''EGGS HEATING - PHASES'''
        # preallocate variables for phase
        self.phase_ch1_turns =          float(0)
        self.phase_phaser_turns_arr =   np.zeros((2, 3), dtype=float)


        '''EGGS HEATING - CONFIG'''
        # convert attenuation from dB to machine units
        self.att_eggs_heating_mu = att_to_mu(self.att_eggs_heating_db * dB)

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


        '''EGGS HEATING - QUANTUM JUMPS'''
        # get mean of sideband frequencies for thresholding
        self.freq_sideband_readout_mean_ftw = np.mean(self.freq_sideband_readout_ftw_list)
        # store sub_repetitions results for real-time monitoring
        self._quantum_jump_monitor_rsb = np.zeros(self.sub_repetitions, dtype=np.int32)
        self._quantum_jump_monitor_bsb = np.zeros(self.sub_repetitions, dtype=np.int32)

        # todo: coherent_c interpolator

        # create config for quantum jumps (just an RSB and a BSB)
        self.config_quantum_jumps = np.zeros((2, np.shape(self.config_eggs_heating_list)[1]))
        self.config_quantum_jumps[:, :] = self.config_eggs_heating_list[0]
        self.config_quantum_jumps[:, 0] = self.freq_sideband_readout_ftw_list
        # create index iterators for quantum jump data structures
        self._results_quantum_jumps_idx = 0
        self._dataset_quantum_jumps_idx = 0

        self.max_num_quantum_jumps = 8  # set max number of possible peaks to prevent boobooing
        self.peak_sb_threshold_min_frac = np.int32(0.10 * self.sub_repetitions) # min threshold for peak detection
        self.peak_sb_threshold_max_frac = np.int32(0.95 * self.sub_repetitions) # max threshold for peak detection


        '''EGGS HEATING - AMPLITUDE CALIBRATION'''
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points = self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve =  Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # TMP REMOVE: MAKE SURE SIDEBAND AMPLITUDES ARE SCALED CORRECTLY FOLLOWING USER INPUT SPECS
        # todo: move to a phaser internal function
        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            for i, (_, carrier_freq_hz, secular_freq_hz, _, _, _, _, _) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz =    (np.array([-secular_freq_hz, secular_freq_hz]) + carrier_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac =        ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac
                # TMP FIX: MAKE SURE SCALED POWER FOLLOWS SPECIFICATIONS OF RSB AND BSB PCT
                # scaled_power_pct =                                          (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                #                                                              ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                scaled_power_pct =  (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                                    ((self.ampl_eggs_heating_rsb_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                # update configs and convert amplitude to frac
                self.config_eggs_heating_list[i, [3, 4, 5]] = np.array([scaled_power_pct[0],
                                                                        scaled_power_pct[1],
                                                                        self.ampl_eggs_dynamical_decoupling_pct]) / 100.


        '''EGGS HEATING - EGGS RF CONFIGURATION'''
        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling:                self.config_eggs_heating_list[:, 5] = 0.

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
        self.time_pulse_shape_rolloff_mu =          self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)

        # todo: document
        self.t_max_phaser_update_rate_mu =          25 * self.phaser_eggs.t_sample_mu
        self.time_pulse_shape_sample_mu =           self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
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
        self.ampl_pulse_shape_frac_list =           np.tile(np.array([self.ampl_eggs_heating_rsb_pct,
                                                                      self.ampl_eggs_heating_bsb_pct,
                                                                      ampl_dd_pct]) / 100.,
                                                            self.num_pulse_shape_samples).reshape(-1, 3)

        # calculate windowing values
        if self.type_pulse_shape == 'sine_squared':
            # calculate sine squared window
            self.ampl_window_frac_list =    np.power(np.sin(
                                                    (np.pi / (2. * self.num_pulse_shape_samples)) *
                                                    np.linspace(1, self.num_pulse_shape_samples, self.num_pulse_shape_samples)),
                                            2)
        elif self.type_pulse_shape == 'error_function':
            raise Exception('Error: error function window not implemented')
        else:
            raise Exception('Error: idk, some window problem')

        # ensure window array has correct dimensions required
        self.ampl_window_frac_list =                np.array([self.ampl_window_frac_list]).transpose()

        # apply window to pulse shape
        self.ampl_pulse_shape_frac_list *=          self.ampl_window_frac_list
        self.ampl_pulse_shape_reverse_frac_list =   self.ampl_pulse_shape_frac_list[::-1]

        # create data structures to hold pulse shaping DMA sequences
        self.phaser_dma_handle_pulseshape_rise = (0, np.int64(0), np.int32(0), False)
        self.phaser_dma_handle_pulseshape_fall = (0, np.int64(0), np.int32(0), False)

    def _prepare_psk(self):
        """
        Calculate and prepare timings for PSK.
        """
        # create config holder for dynamical decoupling PSK; holds time_mu and phase in turns
        self.config_dynamical_decoupling_psk_list = np.zeros((self.num_dynamical_decoupling_phase_shifts + 1, 2), dtype=np.int64)
        self.config_dynamical_decoupling_psk_list[1::2] = 1

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


    # MAIN SEQUENCE
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
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        """
        todo: document
        """
        self.core.break_realtime()

        # get custom sequence handles
        self.phaser_dma_handle_pulseshape_rise = self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        self.phaser_dma_handle_pulseshape_fall = self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()

        # set relevant kernel variables
        _loop_iter = 0  # used to check_termination more frequently


        '''MAIN LOOP'''
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

                # run loop
                counts = self.run_loop(freq_readout_ftw, carrier_freq_hz, sideband_freq_hz, ampl_rsb_frac,
                                       ampl_bsb_frac, ampl_dd_frac, phase_rsb_turns, time_readout_mu)

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

                # death detection
                self.rescue_subsequence.detect_death(counts)
                self.core.break_realtime()

                # check termination more frequently in case reps are low
                if (_loop_iter % 50) == 0:
                    self.check_termination()
                    self.core.break_realtime()
                _loop_iter += 1


                '''STORE QUANTUM JUMPS DATA'''
                # store count data in appropriate holder
                if freq_readout_ftw < self.freq_sideband_readout_mean_ftw:
                    self._quantum_jump_monitor_rsb[_subrep_iter] = counts
                else:
                    self._quantum_jump_monitor_bsb[_subrep_iter] = counts
                self.core.break_realtime()


                '''SUB-REPETITION LOGIC'''
                # finish a single sub_rep set
                if _config_iter % 2 == 1:
                    _subrep_iter += 1

                    # continue sub_repetitions
                    if _subrep_iter < self.sub_repetitions:
                        _config_iter -= 1
                    # move onto next EGGS carrier freq
                    else:
                        # reset relevant loop variables
                        _subrep_iter = 0
                        _config_iter += 1

                        '''QUANTUM JUMP SEARCH'''
                        # move on to quantum jumps we have found a peak
                        if self.process_peak() is True:
                            self.run_quantum_jumps(carrier_freq_hz)
                            self.core.break_realtime()

                # move on to next config in sub_rep set
                else:
                    _config_iter += 1

            '''CLEAN UP A REPETITION'''
            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination at the repetition level
            self.check_termination()
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_loop(self, freq_readout_ftw: TInt32, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat,
                 ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat,
                 phase_rsb_turns: TFloat, time_readout_mu: TInt64) -> TInt32:
        """
        Run main experimental loop.
        Does a full EGGS sequence.
        Arguments:
            ***todo***
        Returns:
                TInt32: the number of PMT counts read out.
        """
        '''CONFIGURE'''
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
        self.phaser_eggs.phaser_setup(self.att_eggs_heating_mu, self.att_eggs_heating_mu)

        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.core_dma.playback_handle(self.phaser_dma_handle_pulseshape_rise)

        # EGGS - RUN
        self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

        # EGGS - STOP
        self.core_dma.playback_handle(self.phaser_dma_handle_pulseshape_fall)
        self.phaser_eggs.phaser_stop()

        '''READOUT'''
        self.sidebandreadout_subsequence.run_time(time_readout_mu)
        self.readout_subsequence.run_dma()
        counts = self.readout_subsequence.fetch_count()

        '''CLEAN UP'''
        # resuscitate ion
        self.rescue_subsequence.resuscitate()
        return counts


    '''
    QUANTUM JUMPS
    '''
    @kernel(flags={"fast-math"})
    def process_peak(self) -> TBool:
        """
        Checks whether a peak has been found using the data from a single sub_rep period.
        A "peak" is considered found if the sideband ratio exceeds the given threshold and the
        constituent RSB and BSB values fall within a threshold range to prevent false positives
        (e.g. from ion death or mass change).
        Returns:
            TBool   : True if peak detected, False otherwise.
        """
        # initialize holder variables
        rsb_total = 0
        bsb_total = 0

        # calculate running total of D-state probabilities
        for i in range(self.sub_repetitions):
            if self._quantum_jump_monitor_rsb[i] <= self.count_threshold:
                rsb_total += 1
            if self._quantum_jump_monitor_bsb[i] <= self.count_threshold:
                bsb_total += 1
        self.core.break_realtime()

        # return peak criteria check
        # note: multiply instead of divide for efficiency and to avoid divide-by-zero errors
        return ((rsb_total > self.peak_sb_threshold_min_frac) and (rsb_total < self.peak_sb_threshold_max_frac) and
                (bsb_total > self.peak_sb_threshold_min_frac) and (bsb_total < self.peak_sb_threshold_max_frac) and
                (rsb_total > bsb_total * self.quantum_jump_threshold_phonon))

    @kernel(flags={"fast-math"})
    def run_quantum_jumps(self, carrier_freq_hz: TFloat) -> TNone:
        """
        Continuously take data at a single carrier frequency to look
        for quantum jumps.
        Arguments:
            carrier_freq_hz (TFloat): the carrier frequency to run quantum jumps at (in Hz).
        """
        # prepare quantum jumps data structures on host-side
        self._prepare_quantum_jumps()

        # run given number of quantum jumps
        for idx_jump in range(self.num_quantum_jumps):

            # loop over configuration list (i.e. RSB, BSB)
            for config_vals in self.config_quantum_jumps:
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                sideband_freq_hz =  config_vals[2]
                ampl_rsb_frac =     config_vals[3]
                ampl_bsb_frac =     config_vals[4]
                ampl_dd_frac =      config_vals[5]
                phase_rsb_turns =   config_vals[6]
                time_readout_mu =   np.int64(config_vals[7])
                self.core.break_realtime()

                # start jumps
                counts = self.run_loop(freq_readout_ftw, carrier_freq_hz, sideband_freq_hz, ampl_rsb_frac,
                                       ampl_bsb_frac, ampl_dd_frac, phase_rsb_turns, time_readout_mu)

                # update results for quantum jumps
                self.update_quantum_jump_results(freq_readout_ftw, counts, carrier_freq_hz)
                self.core.break_realtime()

    @rpc
    def _prepare_quantum_jumps(self) -> TNone:
        """
        Prepare data structures for running quantum jumps.
        """
        # check to see if we're getting too many peaks
        if self._dataset_quantum_jumps_idx >= self.max_num_quantum_jumps:
            raise Exception("Error - too many peaks for quantum jumps: {:d}.".format(self._dataset_quantum_jumps_idx))

        # create new dataset to hold quantum jumps results
        self.set_dataset("results_quantum_jumps_{:d}".format(self._dataset_quantum_jumps_idx),
                         np.zeros((self.num_quantum_jumps * 2, 3)))
        self._dataset_quantum_jumps_idx += 1

        # reset iterators for quantum jumps
        self._results_quantum_jumps_idx = 0

    @rpc(flags={"async"})
    def update_quantum_jump_results(self, *args) -> TNone:
        """
        Store results in the "quantum jumps" dataset.
        Arguments:
            *args: result values to store in dataset.
        """
        # note: use self._dataset_quantum_jumps_idx - 1 since we increment it immediately after creation
        self.mutate_dataset("results_quantum_jumps_{:d}".format(self._dataset_quantum_jumps_idx - 1),
                            self._results_quantum_jumps_idx, np.array(args))
        self._results_quantum_jumps_idx += 1


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

    def analyze_experiment(self):
        # should be backward-compatible with experiment which had no sub-repetitions
        try:
            sub_reps = self.sub_repetitions
        except Exception as e:
            sub_reps = 1

        # create lists for a ch1 sweep
        ch1_turns_sweep_list = []
        phonons_ch1_sweep_list = []

        # determine if a ch1 sweep occurred
        ch1_sweep_bool = len(self.phase_eggs_heating_ch1_turns_list) > 1
        rsb_sweep_bool = len(self.phase_eggs_heating_rsb_turns_list) > 1
        turns_sweep_bool = ch1_sweep_bool or rsb_sweep_bool

        # print results
        print("\tResults - EGGS Heating:")

        # sweep over ch1_turns
        for ch1_turns in self.phase_eggs_heating_ch1_turns_list:

            # handle errors from data processing
            try:
                # convert dataset to array
                dataset_tmp = np.array(self.results)
                dataset = np.reshape(dataset_tmp[np.where(dataset_tmp[:, 5] == ch1_turns), :],
                                     (-1, dataset_tmp.shape[1]))

                ## determine scan type
                col_unique_vals = np.array([len(set(col)) for col in dataset.transpose()])
                # convert unique count to dataset index and order in decreasing value (excluding PMT counts)
                if np.argsort(-col_unique_vals)[0] != 1:
                    sorting_col_num = np.argsort(-col_unique_vals)[0]
                else:
                    sorting_col_num = np.argsort(-col_unique_vals)[1]

                # ensure we actually have a scan, and not some rubbish
                if col_unique_vals[sorting_col_num] <= 1:
                    continue

                ## convert results to sideband ratio
                ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
                                                                                               1, 0,
                                                                                               self.repetitions,
                                                                                               sub_reps)

                # get the phonons and instantiate the fitter class
                phonons = convert_ratios_to_coherent_phonons(ratios)
                fitter = fitSincGeneric()

                # format arguments for applet plotting
                if ch1_sweep_bool:
                    ccb_command = (f'$python -m LAX_exp.applets.plot_matplotlib '
                                   f'temp.plotting.results_eggs_heating_quantum_jumps_{ch1_turns}')
                    group = 'plotting.eggs_heating.ch1_sweep'
                    dataset_name = f'temp.plotting.results_eggs_heating_quantum_jumps_{ch1_turns}'
                    applet_name = f"EGGS Heating - Quantum Jumps - CH1 Turns: {ch1_turns}"
                else:
                    ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_heating_quantum_jumps'
                    group = 'plotting.eggs_heating'
                    dataset_name = f'temp.plotting.results_eggs_heating_quantum_jumps'
                    applet_name = f"EGGS Heating - Quantum Jumps"

                ## process secular or carrier frequency sweep
                if sorting_col_num == 3 or sorting_col_num == 2:

                    fit_x = np.linspace(np.min(scanning_freq_MHz), np.max(scanning_freq_MHz),
                                        10 * len(scanning_freq_MHz))

                    # attempt to fit sinc to data
                    try:
                        # fit swept frequency (secular or carrier) and sidebands
                        fit_params_sweep, fit_err_sweep, _ = fitter.fit(scanning_freq_MHz, phonons)
                        fit_params_rsb, fit_err_rsb, _ = fitter.fit(scanning_freq_MHz, ave_rsb)
                        fit_params_bsb, fit_err_bsb, _ = fitter.fit(scanning_freq_MHz, ave_bsb, -1)

                        # get phonon from fit
                        phonon_n = fit_params_sweep[0]
                        # todo: implement phonon err
                        phonon_err = 0

                        # create arrays for plotting fits
                        fit_y_phonons = fitter.fit_func(fit_x, *fit_params_sweep)
                        fit_y_rsb = fitter.fit_func(fit_x, *fit_params_rsb)
                        fit_y_bsb = fitter.fit_func(fit_x, *fit_params_bsb)
                        # print results to log
                        print("\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_sweep[1] * 1e3,
                                                                          fit_err_sweep[1] * 1e3))

                        # save results to dataset manager for dynamic experiments
                        res_dj = [[phonon_n, phonon_err], [fit_params_sweep, fit_err_sweep]]
                        ch1_turns_sweep_list.append(ch1_turns)
                        phonons_ch1_sweep_list.append(phonon_n)

                        # save results to hdf5 as a dataset
                        if sorting_col_num == 3:
                            self.set_dataset('fit_params_secular', fit_params_sweep)
                            self.set_dataset('fit_err_secular', fit_err_sweep)
                        else:
                            self.set_dataset('fit_params_carrier', fit_params_sweep)
                            self.set_dataset('fit_err_carrier', fit_err_sweep)

                    # if fit fails then ignore plotting of fit by creating array of Nones
                    except Exception as e:
                        fit_y_phonons = [None] * len(fit_x)
                        fit_y_rsb = [None] * len(fit_x)
                        fit_y_bsb = [None] * len(fit_x)
                        res_dj = None

                    # format dictionary of results for plotting with applet
                    ccb_command += ' --num-subplots 3'
                    plotting_results = {'x': [scanning_freq_MHz, scanning_freq_MHz, scanning_freq_MHz],
                                        'y': [ave_rsb, ave_bsb, phonons],
                                        'errors': [std_rsb, std_bsb, [None] * len(std_rsb)],
                                        'fit_x': [fit_x, fit_x, fit_x],
                                        'fit_y': [fit_y_rsb, fit_y_bsb, fit_y_phonons],
                                        'subplot_x_labels': np.array(
                                            ['Frequency (MHz)', 'Frequency (MHz)', 'Frequency (MHz)']),
                                        'subplot_y_labels': np.array(
                                            ['D State Population', 'D State Population', 'Phonons']),
                                        'rid': self.scheduler.rid,
                                        }

                    # determine group to plot applet in
                    if sorting_col_num == 3:
                        group += '.secular'
                    else:
                        group += '.carrier'

                # process sideband readout sweep
                elif sorting_col_num == 0:
                    # todo: implement
                    phonon_err = 0
                    # get the sideband frequencies and prepare for plotting the fit
                    rsb_freqs_MHz, bsb_freqs_MHz, _ = extract_sidebands_freqs(scanning_freq_MHz)
                    fit_x_rsb = np.linspace(np.min(rsb_freqs_MHz), np.max(rsb_freqs_MHz), 1000)
                    fit_x_bsb = np.linspace(np.min(bsb_freqs_MHz), np.max(bsb_freqs_MHz), 1000)

                    try:
                        # try to fit the sidebands
                        fit_params_rsb, fit_err_rsb, fit_rsb = fitter.fit(rsb_freqs_MHz, ave_rsb)
                        fit_params_bsb, fit_err_bsb, fit_bsb = fitter.fit(bsb_freqs_MHz, ave_bsb)
                        fit_y_rsb = fitter.fit_func(fit_x_rsb, *fit_params_rsb)
                        fit_y_bsb = fitter.fit_func(fit_x_bsb, *fit_params_bsb)

                        # get the phonon number and update the list for graphing ch1 sweeps
                        phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
                        ch1_turns_sweep_list.append(ch1_turns)
                        phonons_ch1_sweep_list.append(phonon_n)

                        # save results to hdf5 as a dataset
                        self.set_dataset('fit_params_rsb', fit_params_rsb)
                        self.set_dataset('fit_params_bsb', fit_params_bsb)
                        self.set_dataset('fit_err_rsb', fit_err_rsb)
                        self.set_dataset('fit_err_bsb', fit_err_bsb)

                        # save results to dataset manager for dynamic experiments
                        res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]

                        # print results to log
                        print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2.,
                                                                  float(fit_err_rsb[1]) / 2.))
                        print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2.,
                                                                  float(fit_err_bsb[1]) / 2.))

                    except Exception as e:
                        print("Could not find fit the sidebands")
                        fit_y_rsb = [None] * len(fit_x_rsb)
                        fit_y_bsb = [None] * len(fit_x_bsb)
                        res_dj = None

                    # format dictionary for applet plotting
                    ccb_command += ' --num-subplots 2'
                    plotting_results = {'x': [rsb_freqs_MHz / 2, bsb_freqs_MHz / 2],
                                        'y': [ave_rsb, ave_bsb],
                                        'errors': [std_rsb, std_bsb],
                                        'fit_x': [fit_x_rsb / 2, fit_x_bsb / 2],
                                        'fit_y': [fit_y_rsb, fit_y_bsb],
                                        'subplot_x_labels': np.array(['AOM Frequency (MHz)', 'AOM Frequency (MHz)']),
                                        'subplot_y_labels': np.array(['D State Population', 'D State Population']),
                                        'rid': self.scheduler.rid,
                                        }

                    group += '.sidebands'

                else:
                    # if unknown quantity was scanned ignore analysis
                    plotting_results = {}
                    res_dj = None
                    group = None

                # record values in dataset
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False,
                                 archive=False)
                self.set_dataset(dataset_name, pyon.encode(plotting_results), broadcast=True)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False,
                                 archive=False)

                # create applet
                self.ccb.issue("create_applet", applet_name, ccb_command, group=group)

            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))

        if ch1_sweep_bool:
            try:
                plotting_results = {'x': np.array(ch1_turns_sweep_list)[np.argsort(ch1_turns_sweep_list)],
                                    'y': np.array(phonons_ch1_sweep_list)[np.argsort(ch1_turns_sweep_list)],
                                    'subplot_x_labels': "Channel 1 Phase",
                                    'subplot_y_labels': 'Phonons',
                                    'rid': self.scheduler.rid
                                    }

                self.set_dataset('temp.plotting.results_eggs_heating_quantum_jumps_ch1_sweep',
                                 pyon.encode(plotting_results), broadcast=True)

                self.ccb.issue("create_applet", f"EGGS Heating - Quantum Jumps",
                               '$python -m LAX_exp.applets.plot_matplotlib '
                               'temp.plotting.results_eggs_heating_quantum_jumps_ch1_sweep'
                               ' --num-subplots 1', group=['plotting', 'eggs_heating', 'ch1_sweep'])


            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))

    # def analyze(self):
    #     # should be backward-compatible with experiment which had no sub-repetitions
    #     try:
    #         sub_reps = self.sub_repetitions
    #     except Exception as e:
    #         sub_reps = 1
    #
    #     # handle errors from data processing
    #     try:
    #         # print results
    #         print("\tResults - EGGS Heating:")
    #
    #         # convert dataset to array
    #         dataset = np.array(self.results)
    #
    #         ## determine scan type
    #         # carrier sweep
    #         if len(np.unique(dataset[:, 2])) > 1:       sorting_col_num = 2
    #         # secular frequency
    #         elif len(np.unique(dataset[:, 3])) > 1:     sorting_col_num = 3
    #         else:                                       sorting_col_num = 0
    #
    #         ## convert results to sideband ratio
    #         ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
    #                                                                                        1, 0,
    #                                                                                        self.repetitions, sub_reps)
    #         phonons = convert_ratios_to_coherent_phonons(ratios)
    #
    #         fitter = fitSincGeneric()
    #         rsb_freqs_MHz, bsb_freqs_MHz, _ = extract_sidebands_freqs(scanning_freq_MHz)
    #         fit_params_rsb, fit_err_rsb, fit_rsb = fitter.fit(rsb_freqs_MHz, ave_rsb)
    #         fit_params_bsb, fit_err_bsb, fit_bsb = fitter.fit(bsb_freqs_MHz, ave_bsb)
    #         fit_x_rsb = np.linspace(np.min(rsb_freqs_MHz), np.max(rsb_freqs_MHz), 1000)
    #         fit_x_bsb = np.linspace(np.min(bsb_freqs_MHz), np.max(bsb_freqs_MHz), 1000)
    #         fit_y_rsb = fitter.fit_func(fit_x_rsb, *fit_params_rsb)
    #         fit_y_bsb = fitter.fit_func(fit_x_bsb, *fit_params_bsb)
    #
    #         ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_heating_quantum_jumps'
    #
    #         ## process secular frequency sweep
    #         if sorting_col_num == 3 or sorting_col_num == 2:
    #             fit_params_sweep, fit_err_sweep, _ = fitter.fit(scanning_freq_MHz, phonons)
    #             phonon_n = fit_params_sweep[0]
    #             # todo: implement
    #             phonon_err = 0
    #
    #             # save results to hdf5 as a dataset
    #             if sorting_col_num == 3:
    #                 self.set_dataset('fit_params_secular',  fit_params_sweep)
    #                 self.set_dataset('fit_err_secular',     fit_err_sweep)
    #             else:
    #                 self.set_dataset('fit_params_carrier',  fit_params_sweep)
    #                 self.set_dataset('fit_err_carrier',     fit_err_sweep)
    #
    #
    #             # save results to dataset manager for dynamic experiments
    #             res_dj = [[phonon_n, phonon_err], [fit_params_sweep, fit_err_sweep]]
    #             self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
    #             self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
    #
    #             # print results to log
    #             print("\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_sweep[1] * 1e3, fit_err_sweep[1] * 1e3))
    #
    #             fit_x_phonons = np.linspace(np.min(scanning_freq_MHz), np.max(scanning_freq_MHz), 1000)
    #             fit_y_phonons = fitter.fit_func(fit_x_phonons, *fit_params_sweep)
    #             ccb_command += ' --num-subplots 3'
    #             plotting_results = {'x': [scanning_freq_MHz, scanning_freq_MHz, scanning_freq_MHz],
    #                                 'y': [ave_rsb, ave_bsb, phonons],
    #                                 'errors': [std_rsb, std_bsb, None],
    #                                 'fit_x': [fit_x_rsb, fit_x_bsb, fit_x_phonons],
    #                                 'fit_y': [fit_y_rsb, fit_y_bsb, fit_y_phonons],
    #                                 'subplot_x_labels': ['Frequency (MHz)', 'Frequency (MHz)', 'Frequency (MHz)'],
    #                                 'subplot_y_labels': ['D State Population', 'D State Population', 'Phonons'],
    #                                 'rid': self.scheduler.rid,
    #                                 }
    #
    #
    #         ## process sideband readout sweep
    #         else:
    #             phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
    #             # todo: implement
    #             phonon_err = 0
    #
    #             # save results to hdf5 as a dataset
    #             self.set_dataset('fit_params_rsb',  fit_params_rsb)
    #             self.set_dataset('fit_params_bsb',  fit_params_bsb)
    #             self.set_dataset('fit_err_rsb',     fit_err_rsb)
    #             self.set_dataset('fit_err_bsb',     fit_err_bsb)
    #
    #             # save results to dataset manager for dynamic experiments
    #             res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
    #             self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
    #             self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
    #
    #             # print results to log
    #             print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2., float(fit_err_rsb[1]) / 2.))
    #             print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2., float(fit_err_bsb[1]) / 2.))
    #
    #             ccb_command += ' --num-subplots 2'
    #             plotting_results = {'x': [rsb_freqs_MHz, bsb_freqs_MHz],
    #                                 'y': [ave_rsb, ave_bsb],
    #                                 'errors': [std_rsb, std_bsb],
    #                                 'fit_x': [fit_x_rsb, fit_x_bsb],
    #                                 'fit_y': [fit_y_rsb, fit_y_bsb],
    #                                 'subplot_x_labels': ['Frequency (MHz)', 'Frequency (MHz)'],
    #                                 'subplot_y_labels': ['D State Population', 'D State Population'],
    #                                 'rid': self.scheduler.rid,
    #                                 }
    #
    #         self.set_dataset('temp.plotting.results_eggs_heating_quantum_jumps', pyon.encode(plotting_results),
    #                          broadcast=True)
    #
    #         self.ccb.issue("create_applet", f"EGGS Heating - Quantum Jumps",
    #                        ccb_command)
    #
    #     except Exception as e:
    #         print("Warning: unable to process data.")
    #         print(repr(e))
    #
