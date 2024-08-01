import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)
from sipyco import pyon
from LAX_exp.extensions.physics_constants import *
import matplotlib.pyplot as plt


class EGGSVoltagesCalibration(LAXExperiment, EnvExperiment):
    """
    Calibration: EGGS Voltages

    Examines the resonance for the EGGS feed through and calculates the appropriate amplitude scalings.
    """

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("randomize_config", BooleanValue(default=True))
        self.setattr_argument("sub_repetitions", NumberValue(default=1, ndecimals=0, step=1, min=1, max=500))

        # get subsequences
        self.initialize_subsequence = InitializeQubit(self)
        self.sidebandcool_subsequence = SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence = SidebandReadout(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # readout time - todo: integrate with sideband readout somehow
        self.setattr_argument("time_readout_us", NumberValue(
            default=134.,
            min=1, max=100000, step=1,
            unit="us", scale=1, ndecimals=5
        ), group=self.name)
        # EGGS RF
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list", Scannable(
            default=[
                CenterScan(83.095, 0.8, 0.4, randomize=True),
            ],
            global_min=30, global_max=400, global_step=1,
            unit="MHz", scale=1, ndecimals=6
        ), group='EGGS_Heating.frequencies')

        self.setattr_argument("freq_eggs_heating_secular_one_sb_khz", NumberValue(
            default=1408.00,
            min=0, max=10000, step=1,
            unit="kHz", scale=1, ndecimals=3
        ), group='EGGS_Heating.frequencies')

        self.setattr_argument("freq_eggs_heating_secular_two_sb_khz", NumberValue(
            default=1408.00,
            min=0, max=10000, step=1,
            unit="kHz", scale=1, ndecimals=3
        ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - amplitude
        self.setattr_argument("enable_amplitude_calibration", BooleanValue(default=False),
                              group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",
                              NumberValue(default=20., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",
                              NumberValue(default=20., ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("att_eggs_heating_db", NumberValue(default=10., ndecimals=1, step=0.5, min=0, max=31.5),
                              group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_eggs_heating_ms",
                              NumberValue(default=0.4, ndecimals=5, step=1, min=0.000001, max=10000),
                              group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, ndecimals=3
        ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",
                              NumberValue(default=0., ndecimals=3, step=0.1, min=-1.0, max=1.0),
                              group='EGGS_Heating.waveform.time_phase')

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

        # EGGS RF - dynamical decoupling - configuration
        self.setattr_argument("enable_dynamical_decoupling", BooleanValue(default=True),
                              group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",
                              NumberValue(default=0.01, ndecimals=2, step=10, min=0.0, max=99),
                              group='EGGS_Heating.decoupling')

        # EGGS RF - dynamical decoupling - PSK (Phase-shift Keying)
        self.setattr_argument("enable_dd_phase_shift_keying", BooleanValue(default=False),
                              group='EGGS_Heating.decoupling.psk')
        self.setattr_argument("num_dynamical_decoupling_phase_shifts",
                              NumberValue(default=3, ndecimals=0, step=10, min=1, max=100),
                              group='EGGS_Heating.decoupling.psk')

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

        '''SUBSEQUENCE PARAMETERS'''
        # get readout values
        self.freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.time_readout_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in [self.time_readout_us]])

        '''EGGS HEATING - TIMING'''
        self.time_eggs_heating_mu = self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser sample period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        # todo: move to an internal phaser function
        if self.time_eggs_heating_mu % self.phaser_eggs.t_sample_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_sample_multiples = round(self.time_eggs_heating_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_eggs_heating_mu = np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        '''EGGS HEATING - PHASES'''
        # preallocate variables for phase
        self.phase_ch1_turns = np.float(0)
        self.phase_phaser_turns_arr = np.zeros((2, 3), dtype=float)

        '''EGGS HEATING - CONFIG'''
        # convert build arguments to appropriate values and format as numpy arrays

        self.freq_eggs_carrier_hz_list = np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_two_sb_hz_list = np.array(list([self.freq_eggs_heating_secular_two_sb_khz])) * kHz
        self.freq_eggs_secular_hz_one_sb_list = np.array(list([self.freq_eggs_heating_secular_one_sb_khz])) * kHz
        self.freq_eggs_secular_two_sb_hz = self.freq_eggs_heating_secular_two_sb_khz * kHz
        self.freq_eggs_secular_one_sb_hz = self.freq_eggs_heating_secular_one_sb_khz * kHz
        self.phase_eggs_heating_rsb_turns_list = np.array(list(self.phase_eggs_heating_rsb_turns_list))

        # implement frequency sub-repetitions by "multiplying" the eggs frequency
        self.freq_eggs_carrier_hz_list = np.repeat(self.freq_eggs_carrier_hz_list, self.sub_repetitions)

        # create config data structure with amplitude values
        self.config_eggs_heating_list = np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                  len(self.freq_eggs_carrier_hz_list) *
                                                  len(self.freq_eggs_secular_two_sb_hz_list) *
                                                  len(self.phase_eggs_heating_rsb_turns_list) *
                                                  len(self.time_readout_mu_list),
                                                  8), dtype=float)
        # note: sideband readout frequencies are at the end of the
        # meshgrid to support adjacent_sidebands configuration option
        self.config_eggs_heating_list[:, [1, 2, -2, -1, 0]] = np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list,
                                                                                   self.freq_eggs_secular_two_sb_hz_list,
                                                                                   self.phase_eggs_heating_rsb_turns_list,
                                                                                   self.time_readout_mu_list,
                                                                                   self.freq_sideband_readout_ftw_list),
                                                                       -1).reshape(-1, 5)
        self.config_eggs_heating_list[:, [3, 4, 5]] = np.array([self.ampl_eggs_heating_rsb_pct,
                                                                self.ampl_eggs_heating_bsb_pct,
                                                                self.ampl_eggs_dynamical_decoupling_pct]) / 100.

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
            for i, (_, carrier_freq_hz, secular_freq_hz, _, _, _, _, _) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz = (np.array([-secular_freq_hz, secular_freq_hz]) + carrier_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac = ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac
                # TMP FIX: MAKE SURE SCALED POWER FOLLOWS SPECIFICATIONS OF RSB AND BSB PCT
                # scaled_power_pct =                                          (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                #                                                              ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                scaled_power_pct = (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                                    ((self.ampl_eggs_heating_rsb_pct / 100.) / (
                                            transmitted_power_frac[0] + transmitted_power_frac[1])))
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

        ### PULSE SHAPING - AMPLITUDE WINDOW ###
        # create holder object for pulse amplitudes
        self.ampl_pulse_shape_frac_list = np.tile(np.array([self.ampl_eggs_heating_rsb_pct,
                                                            self.ampl_eggs_heating_bsb_pct,
                                                            self.ampl_eggs_dynamical_decoupling_pct]) / 100.,
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
        self.ampl_pulse_shape_frac_list *= self.ampl_window_frac_list
        self.ampl_pulse_shape_reverse_frac_list = self.ampl_pulse_shape_frac_list[::-1]

    def _prepare_psk(self):
        """
        Calculate and prepare timings for PSK.
        """
        # create config holder for dynamical decoupling PSK; holds time_mu and phase in turns
        self.config_dynamical_decoupling_psk_list = np.zeros((self.num_dynamical_decoupling_phase_shifts + 1, 2),
                                                             dtype=np.int64)
        self.config_dynamical_decoupling_psk_list[1::2] = 1

        # divide total eggs heating time into PSK segments
        self.time_psk_delay_mu = np.int64(
            round(self.time_eggs_heating_mu / (self.num_dynamical_decoupling_phase_shifts + 1)))
        # ensure PSK interval time is a multiple of the phaser sample period
        if self.time_psk_delay_mu % self.phaser_eggs.t_sample_mu:
            # round dynamical decoupling PSK interval to the nearest multiple of phaser sample period
            t_sample_multiples = round(self.time_psk_delay_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_psk_delay_mu = np.int64(t_sample_multiples * self.phaser_eggs.t_sample_mu)

        # update dynamical decoupling config list with verified PSK time
        self.config_dynamical_decoupling_psk_list[:, 0] = self.time_psk_delay_mu
        # ensure that psk rate doesn't exceed the shaping time (t_max_phaser_update_rate_mu; about 25 * t_sample_mu)
        assert self.time_psk_delay_mu >= self.t_max_phaser_update_rate_mu, "Error: num_dynamical_decoupling_phase_shifts too high; PSK update rate exceeds max sustained event rate."

        # set appropriate phaser run method for dynamical decoupling PSK
        if self.enable_dd_phase_shift_keying:
            self.phaser_run = self.phaser_run_psk
        else:
            self.phaser_run = self.phaser_run_nopsk

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list) * 3,
                9)

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

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        self.core.break_realtime()

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_eggs_heating_list:
                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw = np.int32(config_vals[0])
                carrier_freq_hz = config_vals[1]
                ampl_rsb_frac = config_vals[3]
                ampl_bsb_frac = config_vals[4]
                ampl_dd_frac = config_vals[5]
                phase_rsb_turns = config_vals[6]
                time_readout_mu = np.int64(config_vals[7])
                self.core.break_realtime()

                self._run_eggs(carrier_freq_hz, self.freq_eggs_secular_one_sb_hz, phase_rsb_turns, freq_readout_ftw,
                               0., ampl_bsb_frac, ampl_dd_frac, time_readout_mu)

                self._run_eggs(carrier_freq_hz, self.freq_eggs_secular_one_sb_hz, phase_rsb_turns, freq_readout_ftw,
                               ampl_rsb_frac, 0., ampl_dd_frac, time_readout_mu)

                self._run_eggs(carrier_freq_hz, self.freq_eggs_secular_two_sb_hz, phase_rsb_turns, freq_readout_ftw,
                               ampl_rsb_frac, ampl_bsb_frac, 0., time_readout_mu)

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        '''CLEANUP'''
        self.core.break_realtime()
        self.phaser_eggs.reset_oscillators()

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
        carrier_freq_hz = self.config_eggs_heating_list[0, 1]
        sideband_freq_hz = self.config_eggs_heating_list[0, 2]
        phase_rsb_turns = self.config_eggs_heating_list[0, 6]
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
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns)
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        self.phase_ch1_turns = (self.phaser_eggs.phase_inherent_ch1_turns +
                                (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_phaser_turns_arr[0, 0] = phase_rsb_turns
        self.phase_phaser_turns_arr[1, 0] = phase_rsb_turns
        # oscillator 1 (BSB)
        self.phase_phaser_turns_arr[0, 1] = (
                                                    sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        self.phase_phaser_turns_arr[1, 1] = (
                                                    sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
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
    def phaser_stop(self) -> TNone:
        """
        Stop the phaser quickly.
        Set maximum attenuation to prevent output leakage.
        """
        # disable eggs phaser output
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.phaser_eggs.t_sample_mu)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.phaser_eggs.t_sample_mu)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.phaser_eggs.t_sample_mu)

        # switch off EGGS attenuators to prevent leakage
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)

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
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 0],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 0],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 1],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 1],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)

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
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 0],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 0],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 1],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 1],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)

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
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 0],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 0],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 1],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 1],
                                                                          clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 2],
                                                                          clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 2],
                                                                          clr=0)

        # first PSK delay period
        delay_mu(self.config_dynamical_decoupling_psk_list[0][0])

        # conduct PSK on carrier
        for dd_config_vals in self.config_dynamical_decoupling_psk_list[1:]:
            # set oscillator 2 (carrier) with phase shift
            with parallel:
                self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                              phase=self.phase_phaser_turns_arr[
                                                                                        0, 2] + (dd_config_vals[
                                                                                                     1] * 0.5), clr=0)
                self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac,
                                                                              phase=self.phase_phaser_turns_arr[
                                                                                        1, 2] + (dd_config_vals[
                                                                                                     1] * 0.5), clr=0)
                delay_mu(dd_config_vals[0])

    @kernel(flags={"fast-math"})
    def _run_eggs(self, carrier_freq_hz, sideband_freq_hz, phase_rsb_turns, freq_readout_ftw,
                  ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat, time_readout_mu):
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # set default flags
        _squeeze_flag = False
        _rsb_qvsa_flag = False
        _bsb_qvsa_flag = False

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
        # set EGGS attenuators
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()

        # EGGS - RUN
        self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

        # EGGS - STOP
        self.phaser_stop()

        '''READOUT'''
        self.sidebandreadout_subsequence.run_time(time_readout_mu)
        self.readout_subsequence.run_dma()
        counts = self.readout_subsequence.fetch_count()

        if ampl_rsb_frac != 0 and ampl_bsb_frac != 0:
            _squeeze_flag = True
        elif ampl_rsb_frac != 0 and ampl_dd_frac != 0:
            _rsb_qvsa_flag = True
        else:
            _bsb_qvsa_flag = True

        # update dataset
        with parallel:
            self.update_results(
                freq_readout_ftw,
                counts,
                carrier_freq_hz,
                sideband_freq_hz,
                phase_rsb_turns,
                time_readout_mu,
                _rsb_qvsa_flag,
                _bsb_qvsa_flag,
                _squeeze_flag
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

    '''
    ANALYSIS
    '''

    def analyze_experiment(self):

        DATASET_KEY_QUADRUPOLE = "calibration.eggs.scaling_coeffs_quadrupole"
        DATASET_KEY_DIPOLE = "calibration.eggs.scaling_coeffs_dipole"

        try:  # should be backward-compatible with experiment which had no sub-repetitions
            sub_reps = self.sub_repetitions
        except Exception as e:
            sub_reps = 1

        calibrations_eggs_scaling_coeffs_quadrupole = pyon.decode(
            self.get_dataset(DATASET_KEY_QUADRUPOLE, default=pyon.encode({})))
        calibrations_eggs_scaling_coeffs_dipole = pyon.decode(
            self.get_dataset(DATASET_KEY_DIPOLE, default=pyon.encode({})))

        # set up target voltages
        target_V_dipole = 0.015
        target_V_sideband = 8
        # values from experiment
        time_readout_s = self.time_readout_us * 1e-6  # assume single readout time!!!
        secular_freq_two_sb_hz = self.freq_eggs_secular_two_sb_hz   ### assume one secular frequency used
        secular_freq_one_sb_hz = self.freq_eggs_secular_one_sb_hz  ### assume one secular frequency used
        secular_freq_one_sb_ang = 2 * np.pi * secular_freq_one_sb_hz
        secular_freq_two_sb_ang = 2 * np.pi * secular_freq_two_sb_hz
        x01 = np.sqrt(hbar / (2 * mCa * secular_freq_one_sb_ang))
        x02 = np.sqrt(hbar / (2 * mCa * secular_freq_two_sb_ang))


        # grab data for each configuration
        sorting_col_num = 2
        dataset = np.array(self.results)

        rsb_qvsa_dataset = dataset[[bool(val) for val in dataset[:, 6]], :]
        bsb_qvsa_dataset = dataset[[bool(val) for val in dataset[:, 7]], :]
        squeezed_dataset = dataset[[bool(val) for val in dataset[:, 8]], :]

        # extract phonon number
        # tuple with (Ratios, ave_rsb, avd_bsb, std_rsb, std_bsb, scan_freqs)
        rsb_qvsa = extract_ratios(rsb_qvsa_dataset, sorting_col_num,
                                  1, 0, self.repetitions, sub_reps)

        rsb_qvsa_phonons = convert_ratios_to_coherent_phonons(rsb_qvsa[0])

        bsb_qvsa = extract_ratios(bsb_qvsa_dataset, sorting_col_num,
                                  1, 0, self.repetitions, sub_reps)

        bsb_qvsa_phonons = convert_ratios_to_coherent_phonons(bsb_qvsa[0])

        squeezed = extract_ratios(squeezed_dataset, sorting_col_num,
                                  1, 0, self.repetitions, sub_reps)

        squeezed_phonons = convert_ratios_to_squeezed_phonons(squeezed[0])

        assert len(bsb_qvsa_phonons) == len(rsb_qvsa_phonons) == len(squeezed_phonons), "Length Mismath"
        dipole_scaling_coeffs = dict()
        quadrupole_scaling_coeffs = dict()

        Vrs = np.zeros(len(rsb_qvsa_phonons))
        Vbs = np.zeros(len(rsb_qvsa_phonons))
        Vds = np.zeros(len(rsb_qvsa_phonons))
        carrier_freqs_mhz = np.zeros(len(rsb_qvsa_phonons))

        # try converting this from for loop to numpy array operations
        for idx, r_phonon in enumerate(rsb_qvsa_phonons):
            b_phonon = bsb_qvsa_phonons[idx]
            s_phonon = squeezed_phonons[idx]
            carrier_freq_mhz = rsb_qvsa[5][idx]
            carrier_freq_hz = rsb_qvsa[5][idx] * 1e6
            carrier_freq_ang = 2 * np.pi * carrier_freq_hz

            numerator_r = np.sqrt(carrier_freq_ang**2 - secular_freq_two_sb_ang**2) * hbar * np.sqrt(
                np.arcsinh(np.sqrt(s_phonon))) * r_phonon ** (1/4) * r0**2
            denominator_r = qe * np.sqrt(time_readout_s * secular_freq_two_sb_ang) * b_phonon ** (1/4) * x02**2
            Vr = numerator_r / denominator_r

            numerator_b = np.sqrt(
                (carrier_freq_ang**2 - secular_freq_two_sb_ang**2)) * hbar * np.sqrt(
                np.arcsinh(np.sqrt(s_phonon))) * b_phonon**(1/4) * r0**2
            denominator_b = qe * np.sqrt(time_readout_s * secular_freq_two_sb_ang) * r_phonon ** (1 / 4) * x02**2
            Vb = numerator_b / denominator_b

            numerator_d = 2 * np.sqrt(
                ((carrier_freq_ang**2 - secular_freq_one_sb_ang**2)**2)*secular_freq_two_sb_ang * (hbar ** 2) * np.sqrt(
                    r_phonon*b_phonon) * r0 ** 2*x02**4)
            denominator_d = np.sqrt(qe ** 2 * time_readout_s * (carrier_freq_ang**2-secular_freq_two_sb_ang**2) * x01**3
                * np.arcsinh(np.sqrt(s_phonon)))*secular_freq_one_sb_ang
            Vd = numerator_d / denominator_d

            # print(Vr, Vb, Vd)

            Vbs[idx] = Vb
            Vrs[idx] = Vr
            Vds[idx] = Vd
            carrier_freqs_mhz[idx] = carrier_freq_mhz

            # if carrier_freq_mhz in dipole_scaling_coeffs.keys():
            #     dipole_scaling_coeffs[carrier_freq_mhz] = np.mean((target_V_dipole / Vd),
            #                                                       dipole_scaling_coeffs[carrier_freq_mhz])
            # else:
            #     dipole_scaling_coeffs[carrier_freq_mhz] = target_V_dipole / Vd
            #
            # if (carrier_freq_mhz + secular_freq_mhz) in quadrupole_scaling_coeffs.keys():
            #     quadrupole_scaling_coeffs[carrier_freq_mhz + secular_freq_mhz] = np.mean((target_V_sideband / Vb),
            #                                                                              quadrupole_scaling_coeffs[
            #                                                                                  carrier_freq_mhz + secular_freq_mhz])
            # else:
            #     quadrupole_scaling_coeffs[carrier_freq_mhz + secular_freq_mhz] = target_V_sideband / Vb
            #
            # if (carrier_freq_mhz - secular_freq_mhz) in quadrupole_scaling_coeffs.keys():
            #     quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = np.mean((target_V_sideband / Vr),
            #                                                                              quadrupole_scaling_coeffs[
            #                                                                                  carrier_freq_mhz - secular_freq_mhz])
            #     quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = target_V_sideband / Vr
            # else:
            #     quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = target_V_sideband / Vr
        #
        # for key in quadrupole_scaling_coeffs.keys():
        #     if key in calibrations_eggs_scaling_coeffs_quadrupole.keys():
        #         calibrations_eggs_scaling_coeffs_quadrupole[key] *= quadrupole_scaling_coeffs[key]
        #     else:
        #         calibrations_eggs_scaling_coeffs_quadrupole[key] = quadrupole_scaling_coeffs
        #
        # for key in dipole_scaling_coeffs.keys():
        #     if key in calibrations_eggs_scaling_coeffs_dipole.keys():
        #         calibrations_eggs_scaling_coeffs_dipole[key] *= dipole_scaling_coeffs[key]
        #     else:
        #         calibrations_eggs_scaling_coeffs_dipole[key] = dipole_scaling_coeffs

        # self.set_dataset(DATASET_KEY_QUADRUPOLE, pyon.encode(calibrations_eggs_scaling_coeffs_quadrupole))
        # self.set_dataset(DATASET_KEY_DIPOLE, pyon.encode(calibrations_eggs_scaling_coeffs_dipole))

        f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w')
        ax.plot(carrier_freqs_mhz - secular_freq_mhz, Vrs, 'r', label="Red Sideband Voltages")
        ax2.plot(carrier_freqs_mhz + secular_freq_mhz, Vbs, 'b', label="Blue Sideband Voltages")
        ax.set_xlim(np.min(carrier_freqs_mhz - secular_freq_mhz), np.max(carrier_freqs_mhz - secular_freq_mhz))
        ax2.set_xlim(np.min(carrier_freqs_mhz + secular_freq_mhz), np.max(carrier_freqs_mhz + secular_freq_mhz))
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        f.legend()
        plt.figure(2)
        plt.plot(carrier_freqs_mhz, Vds, 'black', label="Carrier Voltages")
        plt.legend()
        plt.show()
