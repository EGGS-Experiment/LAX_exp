import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class EGGSHeatingDipoleTest(SidebandCooling.SidebandCooling):
    """
    Experiment: EGGS Heating Dipole Test

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'EGGS Heating Dipole Test'


    def build_experiment(self):
        # EGGS RF scan configuration
        self.setattr_argument("randomize_config",                           BooleanValue(default=True), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([82]),
                                                                                    # ExplicitScan([82, 83, 86.7]),
                                                                                    # CenterScan(85.1, 0.002, 0.001, randomize=True)
                                                                                ],
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([10000]),
                                                                                    # CenterScan(1202, 2, 0.2, randomize=True)
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating')

        # EGGS RF pulse configuration
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating')
        self.setattr_argument("ampl_eggs_heating_pct",                      NumberValue(default=80, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating')
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=1, ndecimals=5, step=1, min=0.000001, max=10000), group='EGGS_Heating')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=10, ndecimals=1, step=0.5, min=0, max=31.5), group='EGGS_Heating')

        # EGGS RF error correction configuration (e.g. dynamical decoupling)
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=True), group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=20, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.decoupling')

        # get relevant devices
        self.setattr_device('phaser_eggs')
        self.setattr_device('ttl13')
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul0_ch2')
        self.setattr_device('urukul0_ch3')

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # ensure phaser amplitudes sum to less than 100%
        total_phaser_channel_amplitude =                                    self.ampl_eggs_heating_pct + self.ampl_eggs_dynamical_decoupling_pct
        assert total_phaser_channel_amplitude <= 100.,                      "Error: total phaser amplitude exceeds 100%."

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

        ### EGGS HEATING - CONFIG ###
        # create config data structure with amplitude values
        # note: 5 values are [carrier_freq_hz, sideband_freq_hz, rsb_ampl_frac, bsb_ampl_frac, carrier_ampl_frac]
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_eggs_carrier_hz_list) * len(self.freq_eggs_secular_hz_list), 5), dtype=float)
        self.config_eggs_heating_list[:, :2] =                              np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list, self.freq_eggs_secular_hz_list), -1).reshape(-1, 2)
        self.config_eggs_heating_list[:, 2:] =                              np.array([0.4999 * self.ampl_eggs_heating_pct,
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


        # tmp remove
        self.dds_asf = self.urukul0_ch2.amplitude_to_asf(0.35)
        # tmp remove

        # tmp remove: create empty phase holder
        self.dds_pow = np.int32(0)
        self.dds_phase_ch3_turns = 0.
        self.dds_delay_ch3_ns = 2.5
        # tmp remove

        # tmp remove: calculate attenuation register for urukul
        att_val = self.urukul0_cpld.att_to_mu(self.att_eggs_heating_db * dB)
        self.dds_att_reg = np.int32(0x0000FFFF)
        self.dds_att_reg |= ((att_val << (3 * 8)) | (att_val << (2 * 8)))
        # tmp remove

        # run preparations for sideband cooling
        super().prepare_experiment()

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

        # tmp remove: set urukul phase mode
        # with parallel:
        self.urukul0_ch2.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.ttl13.off()
        self.core.break_realtime()
        # tmp remove


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

        # tmp remove: dipole heating sequence
        with self.core_dma.record('_DDS_DIPOLE'):
            # set dds attenuation here - ensures that dds channel will have correct attenuation
            self.urukul0_cpld.set_all_att_mu(self.dds_att_reg)
            self.urukul0_cpld.set_profile(0)

            # reset signal phase
            # with parallel:
            self.urukul0_ch2.set_cfr1(phase_autoclear=1)
            self.urukul0_ch3.set_cfr1(phase_autoclear=1)
            # latch phase reset
            self.urukul0_cpld.io_update.pulse_mu(8)

            # tickle for given time
            with parallel:
                self.urukul0_cpld.cfg_switches(0b1100)
                self.ttl13.on()
            delay_mu(self.time_eggs_heating_mu)
            self.urukul0_cpld.cfg_switches(0b0000)


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        _handle_dds_dipole = self.core_dma.get_handle('_DDS_DIPOLE')
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

                # tmp remove: set up urukul channel frequencies and phase
                freq_dds_ftw = self.urukul0_ch2.frequency_to_ftw(sideband_freq_hz)
                self.dds_pow = self.urukul0_ch2.turns_to_pow((self.dds_delay_ch3_ns * ns * sideband_freq_hz) + self.dds_phase_ch3_turns)
                # tmp remove

                # tmp remove: set urukul0_ch2 and urukul0_ch3 frequencies
                # with parallel:
                self.urukul0_ch2.set_mu(freq_dds_ftw, pow_=0, asf=self.dds_asf, phase_mode=PHASE_MODE_ABSOLUTE, profile=0)
                self.urukul0_ch3.set_mu(freq_dds_ftw, pow_=self.dds_pow, asf=self.dds_asf, phase_mode=PHASE_MODE_ABSOLUTE, profile=0)
                # tmp remove: set urukul0_ch2 and urukul0_ch3 frequencies


                # sweep 729nm readout frequency
                for freq_ftw in self.freq_readout_ftw_list:

                    # set frequency
                    self.qubit.set_mu(freq_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                    self.core.break_realtime()

                    # initialize ion in S-1/2 state
                    self.initialize_subsequence.run_dma()

                    # sideband cool
                    self.sidebandcool_subsequence.run_dma()

                    # run eggs heating (dds dipole, calibration test)
                    self.core_dma.playback_handle(_handle_dds_dipole)

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

        # tmp remove
        self.ttl13.off()
        # tmp remove
