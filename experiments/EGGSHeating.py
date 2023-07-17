import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class HeatingRate(SidebandCooling.SidebandCooling):
    """
    Experiment: EGGS Heating

    Measures the heating rate by doing sideband cooling, then waiting
    a variable amount of time before readout.
    """
    name = 'EGGS Heating'


    def build_experiment(self):
        # EGGS RF scan configuration
        self.setattr_argument("freq_eggs_heating_mhz_carrier_list",         Scannable(
                                                                                default=CenterScan(85.1, 0.020, 0.001, randomize=True),
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ), group='EGGS_Heating')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                # default=CenterScan(1410, 0.001, 0.0001, randomize=True),
                                                                                default=ExplicitScan([1716]),
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating')

        # EGGS RF pulse configuration
        self.setattr_argument("enable_amplitude_calibration",               BooleanValue(default=False), group='EGGS_Heating')
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=2, ndecimals=5, step=1, min=0.000001, max=10000), group='EGGS_Heating')
        self.setattr_argument("ampl_eggs_heating_pct",                      NumberValue(default=90, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating')
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=10, ndecimals=1, step=0.5, min=0, max=31.5), group='EGGS_Heating')

        # EGGS RF error correction configuration (e.g. dynamical decoupling)
        self.setattr_argument("enable_dynamical_decoupling",                BooleanValue(default=False), group='EGGS_Heating.decoupling')
        self.setattr_argument("ampl_eggs_dynamical_decoupling_pct",         NumberValue(default=5, ndecimals=2, step=10, min=0.01, max=99), group='EGGS_Heating.decoupling')

        # get devices
        self.setattr_device('phaser_eggs')

        # run regular sideband cooling build
        super().build_experiment()


    def prepare_experiment(self):
        # todo: check input
        # ensure phaser amplitudes sum to less than 100%
        # min_time_length =                                                   len(list(self.time_min_sideband_cooling_us_list))
        # max_time_length =                                                   len(list(self.time_max_sideband_cooling_us_list))
        # modes_length =                                                      len(list(self.freq_sideband_cooling_mhz_list))
        # assert min_time_length == max_time_length == modes_length

        ### EGGS HEATING - TIMING ###
        self.time_eggs_heating_mu =                                         self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser frame period
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        if self.time_eggs_heating_mu % self.phaser_eggs.t_frame_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_frame_multiples = round(self.time_eggs_heating_mu / self.phaser_eggs.t_frame_mu + 0.5)
            self.time_eggs_heating_mu = np.int64(self.phaser_eggs.t_frame_mu * t_frame_multiples)

        ### EGGS HEATING - FREQUENCIES ###
        self.freq_eggs_carrier_hz_list =                                    np.array(list(self.freq_eggs_heating_mhz_carrier_list)) * MHz - (self.phaser_eggs.freq_center_mhz * MHz)
        self.freq_eggs_secular_hz_list =                                    np.array(list(self.freq_eggs_heating_secular_mhz_list)) * MHz

        ### EGGS HEATING - CONFIG/CALIBRATION ###
        # create config data structure
        # note: 5 values are [center_offset_freq_hz, sideband_freq_hz, rsb_ampl_pct, bsb_ampl_pct, carrier_ampl_pct]
        self.config_eggs_heating_list =                                     np.zeros((len(self.freq_eggs_carrier_hz_list) * len(self.freq_eggs_secular_hz_list), 5), dtype=float)
        self.config_eggs_heating_list[:, :2] =                              np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list, self.freq_eggs_secular_hz_list), -1).reshape(-1, 2)

        # interpolate calibration dataset
        # note: choose 1D interpolator since it ensures smoothness at each point
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                                 self.get_dataset('calibration.eggs.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                                  Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # calculate calibrated eggs sidebands amplitudes
        # todo: implement boolean check for enable_amplitude_calibration
        # todo: correctly set amplitudes
        for i, config_freqs in enumerate(self.config_eggs_heating_list):
            # get relevant frequency values
            carrier_freq, secular_freq = config_freqs[:2]
            carrier_freq += self.freq_eggs_heating_center_hz
            # tmp remove
            print(ampl_calib_curve([(carrier_freq - secular_freq) / MHz, (carrier_freq + secular_freq) / MHz]))
            # calibrate amplitude values
            # self.config_eggs_heating_list[i, 2:] = ampl_calib_curve([(carrier_freq - secular_freq) / MHz, (carrier_freq + secular_freq) / MHz])
            # self.config_eggs_heating_list[i, 2:] = [0.499, 0.499]
            self.config_eggs_heating_list[i, 2:] = [0.49, 0.49, self.ampl_eggs_dynamical_decoupling_pct]

        # if dynamical decoupling is disabled, set carrier amplitude to 0.
        if not self.enable_dynamical_decoupling:
            self.config_eggs_heating_list[:, 5] = 0.

        # run preparations for sideband cooling
        super().prepare_experiment()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list) * len(self.freq_readout_ftw_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # run parent initialization
        super().initialize_experiment()
        # self.core.break_realtime()

        # set attenuations for phaser outputs
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        _handle_eggs_heating = self.core_dma.get_handle('_EGGS_HEATING')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep eggs rf configurations
            for config_vals in self.config_eggs_heating_list:

                # extract values from config list
                carrier_offset_freq_hz =    config_vals[0]
                sideband_freq_hz =          config_vals[1]
                ampl_rsb_frac =             config_vals[2]
                ampl_bsb_frac =             config_vals[3]
                ampl_dd_frac =              config_vals[4]
                self.core.break_realtime()

                # configure EGGS tones
                self.phaser_configure(carrier_offset_freq_hz, sideband_freq_hz)

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
                            (carrier_offset_freq_hz + self.freq_eggs_heating_center_hz) / MHz,
                            sideband_freq_hz
                        )
                        self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

        # CLEANUP
        self.core.break_realtime()
        # reset all oscillator frequencies and amplitudes
        self.phaser_eggs.reset_oscillators()


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_offset_freq_hz: TFloat, sideband_freq_hz: TFloat):
        """
        Configure the tones on phaser for EGGS.
        Puts an RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.
        Arguments:
            carrier_offset_freq_hz  (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
        """
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_offset_freq_hz)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_offset_freq_hz)
        # strobe updates for both channels
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.duc_stb()

        # set sideband frequencies
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz)

        # set carrier frequency (i.e. 0 Hz) for dynamical decoupling
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)

    @kernel(flags={"fast-math"})
    def phaser_run(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
        """
        todo: document
        Puts an RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the carrier sideband amplitude (as a decimal fraction).
        """
        # todo: see if we can set oscillators parallel
        # activate eggs heating output - channel 0
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=0., clr=0)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=0., clr=0)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=0., clr=0)

        # activate eggs heating output - channel 1
        # todo: get phase from phaser_eggs
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=0., clr=0)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=0., clr=0)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=0., clr=0)

        # leave eggs heating running
        delay_mu(self.time_eggs_heating_mu)

        # disable eggs heating output - channel 0
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., clr=1)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., clr=1)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., clr=1)

        # activate eggs heating output - channel 1
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., clr=1)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., clr=1)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., clr=1)
