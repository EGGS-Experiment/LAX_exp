import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class EGGSHeating(SidebandCooling.SidebandCooling):
    """
    Experiment: EGGS Heating

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'EGGS Heating'


    def build_experiment(self):
        # EGGS RF scan configuration
        self.setattr_argument("randomize_config",                           BooleanValue(default=True), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([82, 83, 86.7]),
                                                                                    CenterScan(85.1, 0.002, 0.001, randomize=True)
                                                                                ],
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([1200]),
                                                                                    CenterScan(1202, 2, 0.2, randomize=True)
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating')

        # EGGS RF pulse configuration
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating')
        self.setattr_argument("ampl_eggs_heating_pct",                      NumberValue(default=80, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating')
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=1, ndecimals=5, step=1, min=0.000001, max=10000), group='EGGS_Heating')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=20, ndecimals=1, step=0.5, min=0, max=31.5), group='EGGS_Heating')

        # EGGS RF error correction configuration (e.g. dynamical decoupling)
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=True), group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=20, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.decoupling')

        # get relevant devices
        self.setattr_device('phaser_eggs')

        # run regular sideband cooling build
        super().build_experiment()

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
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_eggs_carrier_hz_list) * len(self.freq_eggs_secular_hz_list), 5), dtype=float)
        self.config_eggs_heating_list[:, :2] =                              np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list, self.freq_eggs_secular_hz_list), -1).reshape(-1, 2)
        self.config_eggs_heating_list[:, 2:] =                              np.array([0.4999 * self.ampl_eggs_heating_pct,
                                                                                      0.4999 * self.ampl_eggs_heating_pct,
                                                                                      self.ampl_eggs_dynamical_decoupling_pct]) / 100.
        # todo: integrate config with readout frequency

        ### EGGS HEATING - CALIBRATION ###
        # interpolate calibration dataset
        # note: we choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                                 self.get_dataset('calibration.eggs.transmission.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                                  Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # calculate calibrated eggs sidebands amplitudes
        if self.enable_amplitude_calibration:
            for i, (carrier_freq_hz, secular_freq_hz, _, _, _) in enumerate(self.config_eggs_heating_list):
                # convert frequencies to absolute units in MHz
                rsb_freq_mhz, bsb_freq_mhz = (np.array([-secular_freq_hz, secular_freq_hz]) + carrier_freq_hz) / MHz
                # get normalized transmission through system
                transmitted_power_frac = ampl_calib_curve([rsb_freq_mhz, bsb_freq_mhz])
                # adjust sideband amplitudes to have equal power and normalize to ampl_eggs_heating_frac
                scaled_power_pct = np.array([transmitted_power_frac[1], transmitted_power_frac[0]]) * ((self.ampl_eggs_heating_pct / 100.) / (transmitted_power_frac[0] + transmitted_power_frac[1]))
                # update configs and convert amplitude to frac
                self.config_eggs_heating_list[i, 2:] = np.array([scaled_power_pct[0], scaled_power_pct[1], self.ampl_eggs_dynamical_decoupling_pct]) / 100.

        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling:
            self.config_eggs_heating_list[:, 4] = 0.

        # if randomize_config is enabled, completely randomize the sweep order
        # i.e. random readout and EGGS heating parameters each iteration, instead of sweeping 1D by 1D
        if self.randomize_config:
            np.random.shuffle(self.config_eggs_heating_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list) * len(self.freq_readout_ftw_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # note: bulk of this code is copy and pasted from SidebandCooling
        # since we can't call "super().initialize_experiment" in a kernel function

        ### BEGIN COPY AND PASTE FROM SIDEBANDCOOLING ###
        self.core.break_realtime()

        # record subsequences onto DMA
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

        # set attenuations for phaser outputs
        self.core.break_realtime()
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep eggs rf configurations
            for config_vals in self.config_eggs_heating_list:

                # extract values from config list
                carrier_freq_hz =           config_vals[0]
                sideband_freq_hz =          config_vals[1]
                ampl_rsb_frac =             config_vals[2]
                ampl_bsb_frac =             config_vals[3]
                ampl_dd_frac =              config_vals[4]
                self.core.break_realtime()

                # configure EGGS tones
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz)
                self.core.break_realtime()

                # sweep 729nm readout frequency
                for freq_ftw in self.freq_readout_ftw_list:

                    # set frequency
                    self.qubit.set_mu(freq_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                    self.core.break_realtime()

                    # initialize ion in S-1/2 state
                    self.initialize_subsequence.run_dma()

                    # sideband cool
                    self.sidebandcool_subsequence.run_dma()

                    # run eggs heating
                    self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

                    # custom SBC readout
                    self.core_dma.playback_handle(_handle_sbc_readout)

                    # update dataset
                    with parallel:
                        self.update_results(
                            freq_ftw,
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
        # align DUCs of both channels by clearing their phase accumulators
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_cfg(clr_once=1)
        # strobe update register to latch change (does simultaneously for both DUCs)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.duc_stb()

        # activate eggs heating output
        at_mu(self.phaser_eggs.get_next_frame_mu())
        time_start_mu = now_mu()

        # set oscillator 0
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=0., clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0, clr=0)
        # set oscillator 1
        at_mu(time_start_mu + self.phaser_eggs.t_sample_mu)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1, clr=0)
        # set oscillator 2
        at_mu(time_start_mu + 2 * self.phaser_eggs.t_sample_mu)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch0_osc2, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_ch1_osc2, clr=0)

        # leave eggs heating running
        delay_mu(self.time_eggs_heating_mu)

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

    def analyze(self):
        print("\tconfig:")
        print("\t\t{}".format(self.config_eggs_heating_list))

        print("\tosc0:")
        print("\t\tphase ch1 osc0: {:.3f}\n".format(self.phase_ch1_osc0))

        print("\tosc1:")
        print("\t\tphase ch0 osc1: {:.3f}".format(self.phase_ch0_osc1))
        print("\t\tphase ch1 osc1: {:.3f}\n".format(self.phase_ch1_osc1))

        print("\tosc2:")
        print("\t\tphase ch0 osc2: {:.3f}".format(self.phase_ch0_osc2))
        print("\t\tphase ch1 osc2: {:.3f}\n".format(self.phase_ch1_osc2))
