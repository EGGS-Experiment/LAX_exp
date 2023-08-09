import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling
# todo: generally migrate everything over to machine units


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
        self.setattr_argument("randomize_config",                           BooleanValue(default=True), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    # ExplicitScan([82, 83, 86.7]),
                                                                                    ExplicitScan([82]),
                                                                                    CenterScan(85.1, 0.002, 0.001, randomize=True)
                                                                                ],
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([1138]),
                                                                                    CenterScan(1138, 2, 0.2, randomize=True)
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.frequencies')

        # EGGS RF pulse configuration
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating.config')
        self.setattr_argument("ampl_eggs_heating_pct",                      NumberValue(default=80, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.config')
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=0.025, ndecimals=5, step=1, min=0.000001, max=10000), group='EGGS_Heating.config')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=10, ndecimals=1, step=0.5, min=3, max=31.5), group='EGGS_Heating.config')

        # EGGS RF - dynamical decoupling
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=True), group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=20, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.decoupling')

        # EGGS RF - pulse shaping
        self.setattr_argument("enable_pulse_shaping",                       BooleanValue(default=True), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",                           EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",                NumberValue(default=100, ndecimals=1, step=100, min=10, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",                NumberValue(default=1500, ndecimals=0, step=100, min=100, max=2000), group='EGGS_Heating.pulse_shaping')

        # get relevant devices
        self.setattr_device('phaser_eggs')

        # tmp remove
        self.setattr_device('ttl13')
        self.setattr_device('ttl16')
        self.setattr_device('ttl17')
        # tmp remove

    def prepare_experiment(self):
        # ensure phaser amplitudes sum to less than 100%
        total_phaser_channel_amplitude =                                    self.ampl_eggs_heating_pct + self.ampl_eggs_dynamical_decoupling_pct
        assert total_phaser_channel_amplitude <= 100.,                      "Error: total phaser amplitude exceeds 100%."

        # run preparations for sideband cooling
        super().prepare_experiment()


        ### EGGS HEATING - TIMING ###
        self.time_eggs_heating_mu =                                         self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser frame period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        # todo: see if it's ok to use sample period instead of frame period
        if self.time_eggs_heating_mu % self.phaser_eggs.t_frame_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_frame_multiples = round(self.time_eggs_heating_mu / self.phaser_eggs.t_frame_mu + 0.5)
            self.time_eggs_heating_mu = np.int64(self.phaser_eggs.t_frame_mu * t_frame_multiples)

        ### EGGS HEATING - FREQUENCIES ###
        self.freq_eggs_carrier_hz_list =                                    np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =                                    np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz

        ### EGGS HEATING - PHASES/TIMING ###
        # preallocate variables for phase
        self.phase_ch0_osc1 = np.float(0)
        self.phase_ch0_osc2 = np.float(0)
        self.phase_ch1_osc0 = np.float(0)
        self.phase_ch1_osc1 = np.float(0)
        self.phase_ch1_osc2 = np.float(0)

        ### EGGS HEATING - CONFIG ###
        # create config data structure with amplitude values
        # note: 5 values are [carrier_freq_hz, sideband_freq_hz, rsb_ampl_frac, bsb_ampl_frac, carrier_ampl_frac]
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_readout_ftw_list) * len(self.freq_eggs_carrier_hz_list) * len(self.freq_eggs_secular_hz_list), 6), dtype=float)
        self.config_eggs_heating_list[:, :3] =                              np.stack(np.meshgrid(self.freq_readout_ftw_list, self.freq_eggs_carrier_hz_list, self.freq_eggs_secular_hz_list), -1).reshape(-1, 3)
        self.config_eggs_heating_list[:, 3:] =                              np.array([0.4999 * self.ampl_eggs_heating_pct,
                                                                                      0.4999 * self.ampl_eggs_heating_pct,
                                                                                      self.ampl_eggs_dynamical_decoupling_pct]) / 100.

        ### EGGS HEATING - CALIBRATION ###
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                                 self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                                  Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            for i, (_, carrier_freq_hz, secular_freq_hz, _, _, _) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz = (np.array([-secular_freq_hz, secular_freq_hz]) + carrier_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac = ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac
                scaled_power_pct = np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) * ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1]))
                # update configs and convert amplitude to frac
                self.config_eggs_heating_list[i, 3:] = np.array([scaled_power_pct[0], scaled_power_pct[1], self.ampl_eggs_dynamical_decoupling_pct]) / 100.

        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling:
            self.config_eggs_heating_list[:, 5] = 0.

        # if randomize_config is enabled, completely randomize the sweep order
        # i.e. random readout and EGGS heating parameters each iteration, instead of sweeping 1D by 1D
        if self.randomize_config:
            np.random.shuffle(self.config_eggs_heating_list)

        # todo: better documentation/implementation
        # configure pulse shaping
        # note: instead of having to deal with adjusting shape, etc., will just add the pulse shaping in addition to the actual pulse
        self._prepare_pulseshape()
        # todo: configure no-op stuff

    def _prepare_pulseshape(self):
        """
        todo: document
        :return:
        """
        ### PULSE SHAPING - TIMING ###
        # convert build variables to units of choice
        self.time_pulse_shape_rolloff_mu =                                  self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)

        # todo: document
        self.t_max_pulse_shape_sample_mu =                                  25 * self.phaser_eggs.t_sample_mu
        self.time_pulse_shape_sample_mu =                                   self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
        # todo: document again lmao
        # todo: add error handling and printouts if there's some sampling problem

        # ensure pulse shaping time is a multiple of the max pulse shape sample rate
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators) is (conservatively) about 1.5 MSPS (i.e. 25 sample periods))
        if self.time_pulse_shape_sample_mu % self.phaser_eggs.t_sample_mu:
            # round pulse shaping sample time up to the nearest multiple of phaser sample period
            t_sample_multiples =                                            round((self.time_pulse_shape_sample_mu / self.phaser_eggs.t_sample_mu) + 0.5)
            self.time_pulse_shape_sample_mu =                               np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        # delay time between successive updates to the pulse envelope, accounts for 2x t_sample_mu delay from having to set 3 oscillators
        self.time_pulse_shape_delay_mu =                                    np.int64(self.time_pulse_shape_sample_mu - 2 * self.phaser_eggs.t_sample_mu)
        # note: calculation of the number of samples accounts for the delay from setting multiple oscillators
        self.num_pulse_shape_samples =                                      np.int32(self.time_pulse_shape_rolloff_mu / (self.time_pulse_shape_sample_mu))


        ### PULSE SHAPING - AMPLITUDE WINDOW ###
        # create holder object for pulse amplitudes
        self.ampl_pulse_shape_frac_list =                                   np.tile(np.array([0.4999 * self.ampl_eggs_heating_pct,
                                                                                              0.4999 * self.ampl_eggs_heating_pct,
                                                                                              self.ampl_eggs_dynamical_decoupling_pct]) / 100., self.num_pulse_shape_samples).reshape(-1, 3)

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
        print('\n\tps sample freq:\t\t{:f} kHz'.format(self.freq_pulse_shape_sample_khz))
        print('\tps sample time (raw):\t\t{:f} ns'.format(self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))))
        print('\tps sample time (aligned):\t{:f} ns'.format(self.time_pulse_shape_sample_mu))
        print('\tps delay time:\t\t{:f} ns'.format(self.time_pulse_shape_delay_mu))
        print('\tps num samples:\t\t{:d}\n'.format(self.num_pulse_shape_samples))

        # print('\n\tampl ps: {}\n\n'.format(self.ampl_pulse_shape_frac_list))
        # raise Exception('stop here idk')

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # note: bulk of this code is copy and pasted from SidebandCooling
        # since we can't call "super().initialize_experiment" in a kernel function

        ### BEGIN COPY AND PASTE FROM SIDEBANDCOOLING ###
        self.core.break_realtime()

        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()

        # record custom readout sequence
        # note: this is necessary since DMA sequences will preserve urukul attenuation register
        with self.core_dma.record('_SBC_READOUT'):
            # set readout waveform for qubit
            self.qubit.set_profile(0)
            self.qubit.set_att_mu(self.att_readout_mu)

            # transfer population to D-5/2 state
            self.rabiflop_subsequence.run()

            # read out fluorescence
            self.readout_subsequence.run()
        ### END COPY AND PASTE FROM SIDEBANDCOOLING ###

        ### PHASER INITIALIZATION ###
        # todo: document better
        # get start
        carrier_freq_hz = self.config_eggs_heating_list[0, 1]
        sideband_freq_hz = self.config_eggs_heating_list[0, 2]
        self.core.break_realtime()
        # configure EGGS tones and set readout frequency
        self.phaser_configure(carrier_freq_hz, sideband_freq_hz)
        self.core.break_realtime()

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
                delay_mu(self.phaser_eggs.t_sample_mu)

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
                delay_mu(self.phaser_eggs.t_sample_mu)

        # set attenuations for phaser outputs
        self.core.break_realtime()
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)
        self.core.break_realtime()

        # tmp remove
        self.ttl13.off()
        self.ttl16.off()
        self.ttl17.off()
        # tmp remove


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_sbc_readout =               self.core_dma.get_handle('_SBC_READOUT')
        _handle_eggs_pulseshape_rise =      self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        _handle_eggs_pulseshape_fall =      self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep eggs rf configurations
            for config_vals in self.config_eggs_heating_list:

                # extract values from config list
                freq_readout_ftw =          np.int32(config_vals[0])
                carrier_freq_hz =           config_vals[1]
                sideband_freq_hz =          config_vals[2]
                ampl_rsb_frac =             config_vals[3]
                ampl_bsb_frac =             config_vals[4]
                ampl_dd_frac =              config_vals[5]
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # run eggs heating
                self.phaser_eggs.reset_duc_phase()
                # todo: hide it away in a method
                self.ttl16.on()
                self.core_dma.playback_handle(_handle_eggs_pulseshape_rise)
                self.ttl17.on()
                self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)
                self.core_dma.playback_handle(_handle_eggs_pulseshape_fall)
                self.phaser_stop()

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        self.readout_subsequence.fetch_count(),
                        carrier_freq_hz,
                        sideband_freq_hz
                    )
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

        # CLEANUP
        self.core.break_realtime()
        # reset all oscillator frequencies and amplitudes
        self.phaser_eggs.reset_oscillators()
        # set max attenuations for phaser outputs to reduce effect of internal noise
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)

        # tmp remove
        self.ttl13.off()
        self.ttl16.off()
        self.ttl17.off()


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat):
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
        """
        # calculate
        # channel 0
        self.phase_ch0_osc1 =           sideband_freq_hz * (self.phaser_eggs.t_sample_mu * ns)
        self.phase_ch0_osc2 =           0.
        # channel 1
        self.phase_ch1_osc0 =           self.phaser_eggs.phase_inherent_ch1_turns +\
                                        ((carrier_freq_hz - sideband_freq_hz) * (self.phaser_eggs.time_latency_ch1_system_ns * ns))
        self.phase_ch1_osc1 =           self.phaser_eggs.phase_inherent_ch1_turns +\
                                        (sideband_freq_hz * (self.phaser_eggs.t_sample_mu * ns)) +\
                                        (carrier_freq_hz + sideband_freq_hz) * (self.phaser_eggs.time_latency_ch1_system_ns * ns)
        # note: extra 0.5 here is to put carrier in dipole config
        self.phase_ch1_osc2 =           self.phaser_eggs.phase_inherent_ch1_turns +\
                                        ((carrier_freq_hz) * (self.phaser_eggs.time_latency_ch1_system_ns * ns)) +\
                                        0.5

        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set sideband frequencies
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # RSB
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # BSB
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # carrier for dynamical decoupling
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_run(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same RSB, BSB, and dynamical decoupling amplitudes for both channels.
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        # activate eggs heating output
        at_mu(self.phaser_eggs.get_next_frame_mu())

        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=0., clr=0)
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

        # tmp remove - commented out for pulse shaping stuff idk
        # disable eggs phaser output
        # with parallel:
        #     self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     delay_mu(self.phaser_eggs.t_sample_mu)
        # with parallel:
        #     self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     delay_mu(self.phaser_eggs.t_sample_mu)
        # with parallel:
        #     self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
        #     self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            # delay_mu(self.phaser_eggs.t_sample_mu)

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
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=0., clr=0)
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

    def analyze(self):
        pass
        # print("\tconfig:")
        # print("\t\t{}".format(self.config_eggs_heating_list))
        #
        # print("\tosc0:")
        # print("\t\tphase ch1 osc0: {:.3f}\n".format(self.phase_ch1_osc0))
        #
        # print("\tosc1:")
        # print("\t\tphase ch0 osc1: {:.3f}".format(self.phase_ch0_osc1))
        # print("\t\tphase ch1 osc1: {:.3f}\n".format(self.phase_ch1_osc1))
        #
        # print("\tosc2:")
        # print("\t\tphase ch0 osc2: {:.3f}".format(self.phase_ch0_osc2))
        # print("\t\tphase ch1 osc2: {:.3f}\n".format(self.phase_ch1_osc2))
