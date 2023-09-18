import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import Squeeze
import LAX_exp.experiments.eggs_heating.EGGSHeating as EGGSHeating

from math import gcd

# todo: generally migrate everything over to machine units
# todo: program phase/profile onto urukul
# todo: synchronize timing for carrier/urukul
# todo: ensure attenuations/power/freq are correctly set


class IonSpectrumAnalyzer(EGGSHeating.EGGSHeating):
    """
    Experiment: Ion Spectrum Analyzer

    ***todo: redocument***
    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'Ion Spectrum Analyzer'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # Ion Spectrum Analyzer
        self.setattr_argument("freq_ionSpecAnal_sideband_offset_khz_list",  Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0]),
                                                                                    CenterScan(0, 5, 0.5, randomize=True)
                                                                                ],
                                                                                global_min=-8000, global_max=8000, global_step=10,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='IonSpectrumAnalyzer')

        # todo: enable squeezing toggling
        self.setattr_argument("enable_squeezing",                       BooleanValue(default=True), group='squeeze')
        self.squeeze_subsequence =                                          Squeeze(self)

    def prepare_experiment(self):
        # ensure phaser amplitudes sum to less than 100%
        total_phaser_channel_amplitude =                                    (self.ampl_eggs_heating_rsb_pct +
                                                                             self.ampl_eggs_heating_bsb_pct +
                                                                             self.ampl_eggs_dynamical_decoupling_pct)
        assert total_phaser_channel_amplitude <= 100.,                      "Error: total phaser amplitude exceeds 100%."

        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list =                                   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list


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
        # convert build arguments to Hz
        self.freq_eggs_carrier_hz_list =                                    np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =                                    np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz
        self.freq_ionSpecAnal_sideband_offset_hz_list =                     np.array(list(self.freq_ionSpecAnal_sideband_offset_khz_list)) * kHz

        # create config data structure with amplitude values
        # note: 6 values are [carrier_freq_hz, sideband_freq_hz, rsb_ampl_frac, bsb_ampl_frac, carrier_ampl_frac, sideband_offset_freq_hz]
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                                                      len(self.freq_eggs_carrier_hz_list) *
                                                                                      len(self.freq_eggs_secular_hz_list) *
                                                                                      len(self.freq_ionSpecAnal_sideband_offset_hz_list),
                                                                                      7), dtype=float)
        self.config_eggs_heating_list[:, [0, 1, 2, -1]] =                   np.stack(np.meshgrid(self.freq_sideband_readout_ftw_list,
                                                                                                 self.freq_eggs_carrier_hz_list,
                                                                                                 self.freq_eggs_secular_hz_list,
                                                                                                 self.freq_ionSpecAnal_sideband_offset_hz_list),
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

        # TMP REMOVE: MAKE SURE SIDEBAND AMPLITUDES ARE SCALED CORRECTLY FOLLOWING USER INPUT SPECS
        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            for i, (_, carrier_freq_hz, secular_freq_hz, _, _, _, offset_freq_hz) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz =                                (np.array([-secular_freq_hz, secular_freq_hz])
                                                                             + carrier_freq_hz + offset_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac =                                    ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac

                # TMP REMOVE
                # TMP FIX: MAKE SURE SCALED POWER FOLLOWS SPECIFICATIONS OF RSB AND BSB PCT
                # scaled_power_pct =                                          (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                #                                                              ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                scaled_power_pct =                                          (np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) *
                                                                            ((self.ampl_eggs_heating_rsb_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1])))
                # update configs and convert amplitude to frac
                self.config_eggs_heating_list[i, [3, 4, 5]] =               np.array([scaled_power_pct[0],
                                                                                      scaled_power_pct[1],
                                                                                      self.ampl_eggs_dynamical_decoupling_pct]) / 100.


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

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_eggs_pulseshape_rise =      self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        _handle_eggs_pulseshape_fall =      self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        # _handle_squeeze =                   self.core_dma.get_handle('_SQUEEZE')
        # _handle_antisqueeze =               self.core_dma.get_handle('_ANTISQUEEZE')
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
                offset_freq_hz =            config_vals[6]
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_psk_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz)
                self.core.break_realtime()
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()


                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()
                # squeeze ion
                # self.core_dma.playback_handle(_handle_squeeze)
                self.squeeze_subsequence.squeeze()


                '''ION SPECTRUM ANALYZER'''
                # PHASER - START/SETUP
                # todo: hide it all away in a method
                self.phaser_eggs.reset_duc_phase()
                self.core_dma.playback_handle(_handle_eggs_pulseshape_rise)

                # PHASER - RUN
                with parallel:
                    # # set TTL for synchronization trigger on a scope
                    # with sequential:
                    #     # delay_mu(self.phaser_eggs.t_output_delay_mu)
                    #     delay_mu(1860)
                    #     self.ttl9.on()
                    # # tmp remove
                    self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

                # PHASER - STOP
                self.core_dma.playback_handle(_handle_eggs_pulseshape_fall)
                self.phaser_stop()


                '''READOUT'''
                # self.core_dma.playback_handle(_handle_antisqueeze)
                self.squeeze_subsequence.antisqueeze()
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        self.readout_subsequence.fetch_count(),
                        carrier_freq_hz,
                        sideband_freq_hz,
                        offset_freq_hz
                    )
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        # CLEANUP
        self.core.break_realtime()
        # reset all oscillator frequencies and amplitudes
        self.phaser_eggs.reset_oscillators()
        # set max attenuations for phaser outputs to reduce effect of internal noise
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)


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
        offset_freq_hz =    self.config_eggs_heating_list[0, 3]
        self.core.break_realtime()

        # configure EGGS tones and set readout frequency; also necessary to ensure phase delays are correctly set
        self.phaser_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz)

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
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, offset_freq_hz: TFloat):
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
        """
        # calculate phase delays between CH0 and CH1
        self.phase_ch1_turns = (self.phaser_eggs.phase_inherent_ch1_turns +
                                (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_ch0_osc0 = self.phase_eggs_heating_rsb_turns
        self.phase_ch1_osc0 = self.phase_eggs_heating_rsb_turns

        # oscillator 1 (BSB)
        self.phase_ch0_osc1 = (sideband_freq_hz + offset_freq_hz) * self.phaser_eggs.t_sample_mu * ns + self.phase_eggs_heating_bsb_turns
        self.phase_ch1_osc1 = (sideband_freq_hz + offset_freq_hz) * self.phaser_eggs.t_sample_mu * ns + self.phase_eggs_heating_bsb_turns

        # oscillator 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_ch0_osc2 = 0.
        self.phase_ch1_osc2 = 0.5
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
        self.phaser_eggs.duc_stb()

        # set sideband frequencies
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz + offset_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz + offset_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz + offset_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz + offset_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)


    # HELPER FUNCTIONS - PSK
    # @kernel(flags={"fast-math"})
    def phaser_psk_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, offset_freq_hz: TFloat):
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
