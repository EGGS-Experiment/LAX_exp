import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.SidebandCooling as SidebandCooling

from math import gcd

# todo: generally migrate everything over to machine units
# todo: program phase/profile onto urukul
# todo: synchronize timing for carrier/urukul
# todo: ensure attenuations/power/freq are correctly set


class EGGSHeating(SidebandCooling.SidebandCooling):
    """
    Experiment: EGGS Heating

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'EGGS Heating'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # EGGS RF scan configuration
        self.setattr_argument("randomize_config",                           BooleanValue(default=True), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([82]),
                                                                                    # ExplicitScan([82, 135, 89.1, 71.3]),
                                                                                    CenterScan(85.1, 0.02, 0.001, randomize=True)
                                                                                ],
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    # ExplicitScan([0.01]),
                                                                                    CenterScan(768.8, 3, 0.2, randomize=True),
                                                                                    ExplicitScan([767.2, 319.2, 1582, 3182])
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - amplitude
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",                  NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",                  NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=4., ndecimals=1, step=0.5, min=0, max=31.5), group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=0.25, ndecimals=5, step=1, min=0.000001, max=10000), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 9, randomize=True)
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",               NumberValue(default=0., ndecimals=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",                       BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",                           EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",                NumberValue(default=100, ndecimals=1, step=100, min=10, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",                NumberValue(default=500, ndecimals=0, step=100, min=100, max=2000), group='EGGS_Heating.pulse_shaping')

        # EGGS RF - dynamical decoupling - configuration
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=True), group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=20., ndecimals=2, step=10, min=0.0, max=99), group='EGGS_Heating.decoupling')

        # EGGS RF - dynamical decoupling - PSK (Phase-shift Keying)
        self.setattr_argument("enable_dd_phase_shift_keying",               BooleanValue(default=False), group='EGGS_Heating.decoupling.psk')
        self.setattr_argument("num_dynamical_decoupling_phase_shifts",      NumberValue(default=3, ndecimals=0, step=10, min=1, max=100), group='EGGS_Heating.decoupling.psk')

        # EGGS RF - dynamical decoupling - active cancellation
        self.setattr_argument("enable_dd_active_cancel",                    BooleanValue(default=False), group='EGGS_Heating.decoupling.activecancel')
        self.setattr_argument("ampl_dd_active_cancel_pct",                  NumberValue(default=35, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.decoupling.activecancel')
        self.setattr_argument("att_dd_active_cancel_db",                    NumberValue(default=10, ndecimals=1, step=0.5, min=3, max=31.5), group='EGGS_Heating.decoupling.activecancel')
        self.setattr_argument("phas_dd_active_cancel_turns",                NumberValue(default=0., ndecimals=3, step=0.1, min=-1, max=1), group='EGGS_Heating.decoupling.activecancel')

        # get relevant devices
        self.setattr_device('phaser_eggs')
        self.setattr_device('urukul1_ch2')

        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')
        # tmp remove

    def prepare_experiment(self):
        # print('yzde123')
        # ensure phaser amplitudes sum to less than 100%
        # total_phaser_channel_amplitude =                                    (self.ampl_eggs_heating_rsb_pct +
        #                                                                      self.ampl_eggs_heating_bsb_pct +
        #                                                                      self.ampl_eggs_dynamical_decoupling_pct)
        # assert total_phaser_channel_amplitude <= 100.,                      "Error: total phaser amplitude exceeds 100%."

        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list =                               self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list


        ### EGGS HEATING - TIMING ###
        self.time_eggs_heating_mu =                                         self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser sample period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        if self.time_eggs_heating_mu % self.phaser_eggs.t_sample_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_sample_multiples =                                            round(self.time_eggs_heating_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_eggs_heating_mu =                                     np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)


        ### EGGS HEATING - PHASES ###
        # preallocate variables for phase
        self.phase_ch1_turns = np.float(0)

        self.phase_ch0_osc0 = np.float(0)
        self.phase_ch0_osc1 = np.float(0)
        self.phase_ch0_osc2 = np.float(0)

        self.phase_ch1_osc0 = np.float(0)
        self.phase_ch1_osc1 = np.float(0)
        self.phase_ch1_osc2 = np.float(0)


        ### EGGS HEATING - CONFIG ###
        # convert build arguments to appropriate values and format as numpy arrays
        self.freq_eggs_carrier_hz_list =                                    np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =                                    np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz
        self.phase_eggs_heating_rsb_turns_list =                            np.array(list(self.phase_eggs_heating_rsb_turns_list))

        # create config data structure with amplitude values
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                                                      len(self.freq_eggs_carrier_hz_list) *
                                                                                      len(self.freq_eggs_secular_hz_list) *
                                                                                      len(self.phase_eggs_heating_rsb_turns_list),
                                                                                      7), dtype=float)
        self.config_eggs_heating_list[:, [0, 1, 2, -1]] =                   np.stack(np.meshgrid(self.freq_sideband_readout_ftw_list,
                                                                                                 self.freq_eggs_carrier_hz_list,
                                                                                                 self.freq_eggs_secular_hz_list,
                                                                                                 self.phase_eggs_heating_rsb_turns_list),
                                                                                     -1).reshape(-1, 4)
        self.config_eggs_heating_list[:, [3, 4, 5]] =                       np.array([self.ampl_eggs_heating_rsb_pct,
                                                                                      self.ampl_eggs_heating_bsb_pct,
                                                                                      self.ampl_eggs_dynamical_decoupling_pct]) / 100.

        ### EGGS HEATING - AMPLITUDE CALIBRATION ###
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                                 self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                                  Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])


        ### EGGS HEATING - EGGS RF CONFIGURATION ###
        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling:                            self.config_eggs_heating_list[:, 5] = 0.

        # if randomize_config is enabled, completely randomize the sweep order
        # i.e. random readout and EGGS heating parameters each iteration, instead of sweeping 1D by 1D
        if self.randomize_config:                                           np.random.shuffle(self.config_eggs_heating_list)

        # configure pulse shaping
        # note: instead of having to deal with adjusting shape, etc., will just add the pulse shaping in addition to the actual pulse
        self._prepare_pulseshape()

        # configure phase-shift keying for dynamical decoupling
        self._prepare_psk()

        # configure active cancellation for dynamical decoupling
        self._prepare_activecancel()

    def _prepare_pulseshape(self):
        """
        todo: document
        :return:
        """
        ### PULSE SHAPING - TIMING ###
        # convert build variables to units of choice
        self.time_pulse_shape_rolloff_mu =                                  self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)

        # todo: document
        self.t_max_phaser_update_rate_mu =                                  25 * self.phaser_eggs.t_sample_mu
        self.time_pulse_shape_sample_mu =                                   self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
        # todo: add error handling and printouts if there's some sampling problem

        # ensure pulse shaping time is a multiple of the max sustained phaser update rate
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators) is (conservatively) about 1.5 MSPS (i.e. 25 sample periods))
        if self.time_pulse_shape_sample_mu % self.t_max_phaser_update_rate_mu:
            # round pulse shaping sample time up to the nearest multiple of phaser sample period
            t_sample_multiples =                                            round((self.time_pulse_shape_sample_mu / self.t_max_phaser_update_rate_mu) + 0.5)
            self.time_pulse_shape_sample_mu =                               np.int64(t_sample_multiples * self.t_max_phaser_update_rate_mu)

        # delay time between successive updates to the pulse envelope, accounts for 2x t_sample_mu delay from having to set 3 oscillators
        self.time_pulse_shape_delay_mu =                                    np.int64(self.time_pulse_shape_sample_mu - 2 * self.phaser_eggs.t_sample_mu)
        # note: calculation of the number of samples accounts for the delay from setting multiple oscillators
        self.num_pulse_shape_samples =                                      np.int32(self.time_pulse_shape_rolloff_mu / (self.time_pulse_shape_sample_mu))


        ### PULSE SHAPING - AMPLITUDE WINDOW ###
        # create holder object for pulse amplitudes
        self.ampl_pulse_shape_frac_list =                                   np.tile(np.array([self.ampl_eggs_heating_rsb_pct,
                                                                                              self.ampl_eggs_heating_bsb_pct,
                                                                                              self.ampl_eggs_dynamical_decoupling_pct]) / 100.,
                                                                                    self.num_pulse_shape_samples).reshape(-1, 3)

        # calculate windowing values
        if self.type_pulse_shape == 'sine_squared':
            # calculate sine squared window
            self.ampl_window_frac_list =                                    np.power(np.sin(
                                                                                (np.pi / (2. * self.num_pulse_shape_samples)) *
                                                                                np.linspace(1, self.num_pulse_shape_samples, self.num_pulse_shape_samples)),
                                                                            2)
        elif self.type_pulse_shape == 'error_function':
            raise Exception('Error: error function window not implemented')
        else:
            raise Exception('Error: idk, some window problem')
        # todo: add other spicy windows
        # ensure window array has correct dimensions required
        self.ampl_window_frac_list =                                        np.array([self.ampl_window_frac_list]).transpose()

        # apply window to pulse shape
        # convert values to machine units (0x3FFF is full scale)
        # self.ampl_window_mu_list =                                      np.int32(self.ampl_window_mu_list * max_amplitude)
        self.ampl_pulse_shape_frac_list *=                                  self.ampl_window_frac_list
        self.ampl_pulse_shape_reverse_frac_list =                           self.ampl_pulse_shape_frac_list[::-1]

        # todo: save select pulse shape values to dataset so we can triple check that everything's ok
        # print('\n\tps sample freq:\t\t{:f} kHz'.format(self.freq_pulse_shape_sample_khz))
        # print('\tps sample time (raw):\t\t{:f} ns'.format(self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))))
        # print('\tps sample time (aligned):\t{:f} ns'.format(self.time_pulse_shape_sample_mu))
        # print('\tps delay time:\t\t{:f} ns'.format(self.time_pulse_shape_delay_mu))
        # print('\tps num samples:\t\t{:d}\n'.format(self.num_pulse_shape_samples))

    def _prepare_psk(self):
        """
        todo: document
        """
        # create config holder for dynamical decoupling PSK; holds time_mu and phase in turns
        self.config_dynamical_decoupling_psk_list =                         np.zeros((self.num_dynamical_decoupling_phase_shifts + 1, 2), dtype=np.int64)
        self.config_dynamical_decoupling_psk_list[1::2] =                   1

        # divide total eggs heating time into PSK segments
        self.time_psk_delay_mu =                                            np.int64(round(self.time_eggs_heating_mu / (self.num_dynamical_decoupling_phase_shifts + 1)))
        # ensure PSK interval time is a multiple of the phaser sample period
        if self.time_psk_delay_mu % self.phaser_eggs.t_sample_mu:
            # round dynamical decoupling PSK interval to the nearest multiple of phaser sample period
            t_sample_multiples =                                            round(self.time_psk_delay_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_psk_delay_mu =                                        np.int64(t_sample_multiples * self.phaser_eggs.t_sample_mu)

        # update dynamical decoupling config list with verified PSK time
        self.config_dynamical_decoupling_psk_list[:, 0] =                   self.time_psk_delay_mu
        # ensure that psk rate doesn't exceed the shaping time (t_max_phaser_update_rate_mu; about 25 * t_sample_mu)
        assert self.time_psk_delay_mu >= self.t_max_phaser_update_rate_mu,  "Error: num_dynamical_decoupling_phase_shifts too high; PSK update rate exceeds max sustained event rate."

        # calculate final delay time as total eggs heating time minus all PSK time
        # note: this pushes all rounding/noninteger problems onto the final delay interval
        # self.time_psk_final_delay_mu =                                      self.time_eggs_heating_mu - np.sum(self.config_dynamical_decoupling_psk_list[:, 0])

        # todo: document
        # set scaling factor to increase resolution of detuned carrier period for PSK period calculation
        self.time_psk_scaling_factor =                                      1e11

        # set appropriate phaser run method for dynamical decoupling PSK
        if self.enable_dd_phase_shift_keying:                               self.phaser_run = self.phaser_run_psk
        else:                                                               self.phaser_run = self.phaser_run_nopsk

    def _prepare_activecancel(self):
        """
        todo: document
        """
        # convert build variables to machine units (where possible)
        self.ampl_dd_active_cancel_asf = self.urukul1_ch2.amplitude_to_asf(self.ampl_dd_active_cancel_pct / 100.)
        self.att_dd_active_cancel_mu = self.urukul1_ch2.cpld.att_to_mu(self.att_dd_active_cancel_db * dB)

        # preallocate variables for phase
        self.phase_activecancel_turns = np.float(0)

        # set appropriate run method for dynamical decoupling active cancellation
        if self.enable_dd_active_cancel:
            self.phaser_activecancel_run = self.phaser_activecancel_run_actual
        else:
            self.phaser_activecancel_run = self.phaser_activecancel_run_pass

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # # todo: set up attenuation for cancellation dds
        # self.core.break_realtime()
        # self.urukul1_ch2.set_att_mu(self.att_dd_active_cancel_mu)
        # self.core.break_realtime()
        # # tmp remove

        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        ### PHASER INITIALIZATION ###
        self.phaser_setup()
        self.core.break_realtime()

        # tmp remove
        self.ttl8.off()
        self.ttl9.off()
        # tmp remove


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_eggs_pulseshape_rise =      self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        _handle_eggs_pulseshape_fall =      self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()


        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep eggs rf configurations
            for config_vals in self.config_eggs_heating_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =          np.int32(config_vals[0])
                carrier_freq_hz =           config_vals[1]
                sideband_freq_hz =          config_vals[2]
                ampl_rsb_frac =             config_vals[3]
                ampl_bsb_frac =             config_vals[4]
                ampl_dd_frac =              config_vals[5]
                phase_rsb_turns =           config_vals[6]
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_psk_configure(carrier_freq_hz, sideband_freq_hz)
                self.core.break_realtime()
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
                # todo: hide it all away in a method
                self.phaser_eggs.reset_duc_phase()
                # tmp remove - integrator hold
                # self.ttl10.on()
                # tmp remove - integrator hold
                self.core_dma.playback_handle(_handle_eggs_pulseshape_rise)

                # EGGS - RUN
                with parallel:
                    self.phaser_activecancel_run()
                    self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

                # EGGS - STOP
                self.core_dma.playback_handle(_handle_eggs_pulseshape_fall)
                with parallel:
                    self.phaser_stop()
                    self.phaser_activecancel_stop()
                # tmp remove - integrator hold
                # self.ttl10.off()
                # tmp remove - integrator hold


                '''READOUT'''
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        self.readout_subsequence.fetch_count(),
                        carrier_freq_hz,
                        sideband_freq_hz,
                        phase_rsb_turns
                    )
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        '''CLEANUP'''
        self.core.break_realtime()
        # reset all oscillator frequencies and amplitudes
        self.phaser_eggs.reset_oscillators()
        # set max attenuations for phaser outputs to reduce effect of internal noise
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)


    # ANALYSIS
    def analyze_experiment(self):

        # get relevant data
        _MU_to_MHz = 2 * 2.32830644e-7
        results_tmp = np.array([self.results])
        readout_freqs = np.array([results_tmp[:.0]]) * _MU_to_MHz
        counts = np.array([results_tmp[:, 1]])
        carrier_freqs_hz = np.array([results_tmp[:,2]]) * _MU_to_MHz
        sideband_freqs_hz = np.array([results_tmp[:,3]])* _MU_to_MHz
        probs = np.zeros(len(counts))

        print(np.shape(result_tmp))


        threshold_list = findThresholdScikit(counts)

        for threshold in threshold_list:
            probs[np.where(counts>threshold)] +=1

        normalized_probs = 1.-probs/len(threshold_list)









    # HELPER FUNCTIONS - PHASER
    @kernel(flags={"fast-math"})
    def phaser_setup(self):
        """
        todo: document
        """
        # todo: document better
        # get starting phase values for pulse shaping ### todo: document better
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

        # set attenuations for phaser outputs
        self.core.break_realtime()
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, phase_rsb_turns: TFloat):
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns)
        """
        # calculate phase delays between CH0 and CH1
        self.phase_ch1_turns =          (self.phaser_eggs.phase_inherent_ch1_turns +
                                         (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_ch0_osc0 =           phase_rsb_turns
        self.phase_ch1_osc0 =           phase_rsb_turns

        # oscillator 1 (BSB)
        self.phase_ch0_osc1 =           (sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        self.phase_ch1_osc1 =           (sideband_freq_hz * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns

        # oscillator 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_ch0_osc2 =           0.
        self.phase_ch1_osc2 =           0.5
        self.core.break_realtime()


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


        # set sideband frequencies
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
    def phaser_stop(self):
        """
        tmp remove
        todo document?
        :return:
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
        # todo: set attenuators?

    @kernel(flags={"fast-math"})
    def phaser_pulseshape_point(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
        """
        todo: document
        :param rsb_ampl:
        :param bsb_ampl:
        :param dd_ampl:
        :return:
        """
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch0_osc0, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2, clr=0)


    # HELPER FUNCTIONS - PSK
    # @kernel(flags={"fast-math"})
    def phaser_psk_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat):
        """
        todo: document

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
        """
        # calculate scaled detuning timings
        time_period_scaled =                                                round(self.time_psk_scaling_factor / (carrier_freq_hz - sideband_freq_hz))
        time_sample_scaled =                                                round(40e-9 * self.time_psk_scaling_factor)

        # calculate lowest common multiple of the scaled detuning period and the scaled sample period
        time_lcm_mu =                                                       round((time_period_scaled * time_sample_scaled) /
                                                                                  gcd(time_period_scaled, time_sample_scaled) *
                                                                                  (1e9 / self.time_psk_scaling_factor))
        # divide total eggs heating time into PSK segments
        time_psk_tmp_delay_mu =                                             round(self.time_eggs_heating_mu / (self.num_dynamical_decoupling_phase_shifts + 1))

        # print(time_period_scaled)
        # print(time_sample_scaled)
        # print('\n')
        # print(time_lcm_mu)
        # print(time_psk_tmp_delay_mu)

        # ensure PSK interval time is very close to a multiple of the carrier detuning period
        if time_psk_tmp_delay_mu % time_lcm_mu:
            # round dynamical decoupling PSK interval to the nearest multiple of phaser sample period
            t_period_multiples =                                            round(time_psk_tmp_delay_mu / time_lcm_mu)
            time_psk_tmp_delay_mu =                                         t_period_multiples * time_lcm_mu

        # update dynamical decoupling config list with new PSK time
        self.config_dynamical_decoupling_psk_list[:, 0] =                   np.int64(time_psk_tmp_delay_mu)

    @kernel(flags={"fast-math"})
    def phaser_run_nopsk(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
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
            self.ttl8.on()
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch0_osc0, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2, clr=0)

        # main eggs pulse
        delay_mu(self.time_eggs_heating_mu)
        self.ttl8.off()

    @kernel(flags={"fast-math"})
    def phaser_run_psk(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
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
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch0_osc0, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2, clr=0)

        # first PSK delay period
        delay_mu(self.config_dynamical_decoupling_psk_list[0][0])

        # conduct PSK on carrier
        for dd_config_vals in self.config_dynamical_decoupling_psk_list[1:]:
            # set oscillator 2 (carrier) with phase shift
            with parallel:
                self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2 + (dd_config_vals[1] * 0.5), clr=0)
                self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2 + (dd_config_vals[1] * 0.5), clr=0)
                # self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1 + (dd_config_vals[1] * 0.5), clr=0)
                # self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1 + (dd_config_vals[1] * 0.5), clr=0)
                delay_mu(dd_config_vals[0])


    # HELPER FUNCTIONS - ACTIVE CANCELLATION
    @kernel(flags={"fast-math"})
    def phaser_activecancel_configure(self, sideband_freq_hz: TFloat):
        """
        todo: document
        """
        # todo: calculate phase delays for the active cancellation DDS channel
        # todo: relative urukul-phaser output delay, system delay, phase argument
        # self.phase_activecancel_turns = self.phaser_eggs.phase_inherent_ch1_turns +\
        #                                 (sideband_freq_hz * (self.phaser_eggs.time_latency_ch1_system_ns) * ns)
        self.phase_activecancel_turns = 0.

        # set up cancellation waveform and profile
        freq_activecancel_ftw = self.urukul1_ch2.frequency_to_ftw(sideband_freq_hz)
        phase_activecancel_pow = self.urukul1_ch2.turns_to_pow(self.phase_activecancel_turns)
        self.urukul1_ch2.set_mu(freq_activecancel_ftw, pow_=phase_activecancel_pow, asf=self.ampl_dd_active_cancel_asf, phase_mode=PHASE_MODE_ABSOLUTE, profile=1)

    @kernel(flags={"fast-math"})
    def phaser_activecancel_run_actual(self):
        """
        todo: document
        """
        # reset signal phase
        delay_mu(350)
        self.urukul1_ch2.cpld.set_profile(1)
        self.urukul1_ch2.set_cfr1(phase_autoclear=1)
        self.urukul1_ch2.cpld.io_update.pulse_mu(8)
        self.urukul1_ch2.sw.off()

    @kernel(flags={"fast-math"})
    def phaser_activecancel_run_pass(self):
        """
        todo: document
        """
        pass

    @kernel(flags={"fast-math"})
    def phaser_activecancel_stop(self):
        """
        todo: document
        """
        # tmp do better/implement better
        self.urukul1_ch2.sw.off()
        # self.urukul1_ch2.set(350 * MHz, amplitude=0.01)
        # self.urukul1_ch2.set_att_mu(0xFF)