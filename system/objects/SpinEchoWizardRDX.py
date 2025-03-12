import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXEnvironment
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

import matplotlib.pyplot as plt
# todo: make display waveform support FFTing to show leakage etc.
# todo: kernel invariants, annotate types & returns


class SpinEchoWizardRDX(LAXEnvironment):
    """
    Helper: Spin Echo Wizard RDX

    Design a spin-echo type pulse sequence for phaser based on user input.
    Output format is designed for Phaser Pulse Shaper object.
    Does not interact with the core device (but requires it for value conversion).
    """
    name = 'Spin Echo Wizard RDX'
    # kernel_invariants = {
    #     "_max_waveforms", "t_max_phaser_update_rate_mu",
    #     "_dma_names", "_dma_handles"
    # }

    def build(self):
        """
        Build SpinEchoWizard and create configuration variables.
        """
        # get relevant devices
        # todo: make it somehow independent of core
        self.setattr_device('core')
        self.setattr_device('phaser_eggs')

        # set max hardware sample rate
        # note: max update rate should be a multiple of 5x the sample period
        # such that each oscillator is deterministically updated
        self.t_max_phaser_update_rate_mu =  5 * self.phaser_eggs.t_sample_mu
        self.num_max_phaser_samples =       190

        '''WAVEFORM CONFIGURATION SEQUENCE'''
        self.sequence_blocks = [
            {
                "oscillator_parameters": [
                    [40., -0.0],    # [ampl_osc0, phas_osc0]
                    [40., 0.0],     # [ampl_osc1, phas_osc1]
                    [20., 0.]       # [ampl_osc2, phas_osc2]
                ],
                "config": {
                    "time_us":          5000,
                    "pulse_shaping":    True,
                    "pulse_shaping_config": {
                        "pulse_shape":          "sine_squared",
                        "pulse_shape_rising":   True,
                        "pulse_shape_falling":  True,
                        "sample_rate_khz":      500,
                        "rolloff_time_us":      100
                    }
                }
            }
        ]

    def _calculate_pulseshape(self, num_oscs: TInt32, pulse_shape: TStr,
                              sample_rate_khz: TFloat, time_rolloff_us: TFloat)\
            -> TTuple([TArray(TInt64, 1), TArray(TFloat, 1)]):
        """
        Calculate desired pulseshape per block based on configured options.
        Arguments:
            num_oscs: the number of oscillators used.
            pulse_shape: the desired pulse shape.
            sample_rate_khz: the pulse-shape update sample rate.
            time_rolloff_us: the total pulse shaping time.
        Returns:
            a tuple of (time_update_arr_mu, amplitude_pulse_shape), where
                time_update_arr_mu is the time to hold/delay after an amplitude update, and
                amplitude_pulse_shape is in fractional units (i.e. between [0., 1.]).
        """
        '''calculate valid update interval and multiple of phaser sample rate'''
        update_interval_mu = self.core.seconds_to_mu(1. / (sample_rate_khz * kHz))
        if update_interval_mu < self.t_max_phaser_update_rate_mu:
            raise ValueError("Waveform sample rate too fast: {:d}ns update interval".format(update_interval_mu))
        # note: (num_oscs - 1) accounts for additional t_sample_mu delay from having to set extra oscillators
        phaser_periods_per_sample = round(update_interval_mu / self.phaser_eggs.t_sample_mu) - (num_oscs - 1)
        update_interval_mu = np.int64(phaser_periods_per_sample * self.phaser_eggs.t_sample_mu)

        '''calculate total number of samples and ensure valid'''
        time_rolloff_mu = self.core.seconds_to_mu(time_rolloff_us * us)
        num_samples = round(time_rolloff_mu / update_interval_mu)
        if num_samples > self.num_max_phaser_samples:
            raise ValueError("Too many points in phaser waveform. Reduce sample rate or rollon time.")
        time_rolloff_mu = np.int64(num_samples * update_interval_mu)

        '''calculate pulse shape (i.e. y-vals)'''
        try:
            # get corresponding pulse shape function
            pulse_shape_func = available_pulse_shapes[pulse_shape]
            ampl_pulse_shape_vals = pulse_shape_func(
                np.arange(num_samples, dtype=np.int64) * update_interval_mu,    # array of time-series values
                time_rolloff_mu
            )
        except KeyError:
            raise ValueError('Invalid pulse shape: {}'.format(pulse_shape))

        time_update_arr_mu = np.ones(num_samples, dtype=np.int64) * update_interval_mu
        return time_update_arr_mu, ampl_pulse_shape_vals

    def _floor_phaser_sample_interval(self, time_us: TFloat) -> TInt64:
        """
        Ensure a given interval is a multiple of the phaser sample period.
        Takes the floor.
        """
        time_mu = self.core.seconds_to_mu(time_us * us)
        time_mu -= time_mu % self.phaser_eggs.t_sample_mu
        return time_mu

    def compile_waveform(self):
        """
        Compile configured waveform for use with PhaserPulseShaper.
        """
        '''VERIFY VALID WAVEFORM SEQUENCE'''
        # get shape of oscillator parameters for each block
        try:
            osc_wav_shapes = np.array([
                np.shape(block["oscillator_parameters"])
                for block in self.sequence_blocks
            ])
        except ValueError:
            raise ValueError("Invalid waveform sequence."
                             "Ensure amplitude and phase are specified for all oscillators in each block.")

        # ensure oscillator parameters have uniform shape across each block
        if np.any(osc_wav_shapes[:, 0] != osc_wav_shapes[0, 0]):
            raise ValueError("Invalid waveform sequence."
                             "All blocks must update the same number of oscillators.")
        elif np.any(osc_wav_shapes[:, 1] != 2):
            raise ValueError("Invalid waveform sequence."
                             "Must specify exactly 2 parameters (ampl and phase) for each oscillator update.")

        # ensure total oscillator amplitudes never exceed max (i.e. < 100%)
        # todo

        # get shapes of blocks
        num_blocks =    len(self.sequence_blocks)
        num_oscs =      osc_wav_shapes[0, 0]

        # create temporary holder objects
        _ampl_arrs =    list(np.zeros((num_oscs, 0), dtype=np.float64))
        _phas_arrs =    list(np.zeros((num_oscs, 0), dtype=np.float64))
        _time_arr =     np.zeros(0, dtype=np.int64)

        # tmp remove - quick workaround for display_waveform
        self.num_blocks = num_blocks
        self.num_oscs = num_oscs
        # tmp remove - quick workaround for display_waveform


        '''BEGIN COMPILATION'''
        for idx_block, _block in enumerate(self.sequence_blocks):

            # get block configuration values
            block_vals =    _block["oscillator_parameters"]
            block_config =  _block["config"]
            block_time_us = block_config["time_us"]

            # set spin-echo/zero-update delay if block_vals is an empty list
            if (block_vals is None) or (len(block_vals) == 0):
                # clear ampl and phase for each osc
                for idx_osc in range(num_oscs):
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], 0.)
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], 0.)
                # set delay time
                _time_arr = np.append(_time_arr, self._floor_phaser_sample_interval(block_time_us))


            # otherwise, configure as normal block
            else:
                # get pulse shaping configuration for block
                pulse_shaping_status =  block_config.get("pulse_shaping", False)
                pulse_shaping_config =  block_config.get("pulse_shaping_config", {})
                pulse_shape_rising =    pulse_shaping_config.get("pulse_shape_rising", False)
                pulse_shape_falling =   pulse_shaping_config.get("pulse_shape_falling", False)


                # implement pulse shaping: rising edge
                if pulse_shaping_status and pulse_shape_rising:

                    # get pulse/window shape of block
                    time_tmp, ampl_window_rising = self._calculate_pulseshape(
                        num_oscs,
                        pulse_shaping_config["pulse_shape"],
                        pulse_shaping_config["sample_rate_khz"],
                        pulse_shaping_config["rolloff_time_us"]
                    )

                    # scale each oscillator's amplitude+phase identically
                    for idx_osc, block_osc_vals in enumerate(block_vals):
                        wind_ampl_tmp = ampl_window_rising * block_osc_vals[0]
                        wind_phas_tmp = np.ones(len(ampl_window_rising)) * block_osc_vals[1]
                        _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)
                    _time_arr = np.append(_time_arr, time_tmp)


                # implement main pulse
                for idx_osc, block_osc_vals in enumerate(block_vals):
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], block_osc_vals[0])
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], block_osc_vals[1])
                _time_arr = np.append(_time_arr, self._floor_phaser_sample_interval(block_time_us))


                # implement pulse shaping: falling edge
                if pulse_shaping_status and pulse_shape_falling:

                    # get pulse/window shape of block
                    time_tmp, ampl_window_rising = self._calculate_pulseshape(
                        num_oscs,
                        pulse_shaping_config["pulse_shape"],
                        pulse_shaping_config["sample_rate_khz"],
                        pulse_shaping_config["rolloff_time_us"]
                    )
                    ampl_window_falling = ampl_window_rising[::-1]

                    # scale each oscillator's amplitude+phase identically
                    for idx_osc, block_osc_vals in enumerate(block_vals):
                        wind_ampl_tmp = ampl_window_falling * block_osc_vals[0]
                        wind_phas_tmp = np.ones(len(ampl_window_falling)) * block_osc_vals[1]
                        _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)
                    _time_arr = np.append(_time_arr, time_tmp)

        # store results
        self._ampl_tmp_arr = np.array(_ampl_arrs)
        self._phas_tmp_arr = np.array(_phas_arrs)
        self._time_tmp_arr = np.array(_time_arr)

    def get_waveform(self) -> TTuple([TArray(TFloat, 2),
                                      TArray(TFloat, 2),
                                      TArray(TInt64, 1)]):
        """
        Retrieve the designed waveform.
        todo: annotate return types
        """
        # note: need to convert ampls from pct to frac
        return (self._ampl_tmp_arr.transpose() / 100.,
                self._phas_tmp_arr.transpose(),
                self._time_tmp_arr)

    def display_waveform(self):
        """
        Display waveform phase and amplitude of each oscillator.
        """
        try:
            # get sequence shape
            # num_blocks, num_oscs, _ = np.shape(self.sequence_blocks)

            # tmp remove - quick workaround for display_waveform
            num_blocks = self.num_blocks
            num_oscs = self.num_oscs
            # tmp remove - quick workaround for display_waveform

            # get cumulative times from pulse delays in _time_tmp_arr
            # note: since timings are post-update delays, need to ensure we start at 0
            _x_axis_time_mu = np.cumsum(np.concatenate(([0], self._time_tmp_arr)))[:-1]
            # convert units from mu to us
            _x_axis_time_us = self.core.mu_to_seconds(_x_axis_time_mu) / us

            # create waveform summary figure
            fig, axs = plt.subplots(num_oscs, 2, layout='constrained')
            fig.suptitle('Spin-Echo Wizard\nTotal pulse time: {:.3f} us'.format(_x_axis_time_us[-1]),
                         fontsize=14)
            # create color list to serve truth and beauty
            colors = plt.cm.rainbow(np.linspace(0, 1, num_oscs))

            # plot waveforms for each oscillator
            for idx in range(num_oscs):
                # waveform plot
                axs[idx, 0].set_title('Amplitude: osc_{:d}'.format(idx))
                axs[idx, 0].set_xlabel('Time (us)')
                axs[idx, 0].set_xlim(0., _x_axis_time_us[-1])
                axs[idx, 0].set_ylabel('Amplitude (%)')
                axs[idx, 0].set_ylim(0., 100.)
                axs[idx, 0].plot(_x_axis_time_us, self._ampl_tmp_arr[idx],
                                 label='osc_{:}'.format(idx), c=colors[idx],
                                 drawstyle="steps-post")

                # phase plot
                axs[idx, 1].set_title('Phase: osc_{:d}'.format(idx))
                axs[idx, 1].set_xlabel('Time (us)')
                axs[idx, 1].set_xlim(0., _x_axis_time_us[-1])
                axs[idx, 1].set_ylabel('Phase (turns)')
                axs[idx, 1].set_ylim(-1., 1.)
                axs[idx, 1].plot(_x_axis_time_us, self._phas_tmp_arr[idx],
                                 label='osc_{:}'.format(idx), c=colors[idx],
                                 drawstyle="steps-post")

            # debugging
            # print(self._ampl_tmp_arr)
            # print(self._phas_tmp_arr)
            # print(self._time_tmp_arr)

            # display figure
            # plt.show(block=False)
            plt.show()

        except Exception as e:
            print(repr(e))


# # for testing purposes
# if __name__ == "__main__":
#     wizard = SpinEchoWizardRDX()
#     wizard.build()
#     wizard.sequence_blocks = [
#         { # block 0
#             "oscillator_parameters": [
#                 [20., -0.0],  # [ampl_osc0, phas_osc0]
#                 [20., 0.0],  # [ampl_osc1, phas_osc1]
#                 [10., 0.]  # [ampl_osc2, phas_osc2]
#             ],
#             "config": {
#                 "time_us": 500,
#                 "pulse_shaping": True,
#                 "pulse_shaping_config": {
#                     "pulse_shape": "sine_squared",
#                     "pulse_shape_rising": True,
#                     "pulse_shape_falling": True,
#                     "sample_rate_khz": 500,
#                     "rolloff_time_us": 100
#                 }
#             }
#         },
#         { # block 1
#             "oscillator_parameters": [
#                 [40., -0.0],  # [ampl_osc0, phas_osc0]
#                 [40., 0.0],  # [ampl_osc1, phas_osc1]
#                 [20., 0.]  # [ampl_osc2, phas_osc2]
#             ],
#             "config": {
#                 "time_us": 250,
#                 "pulse_shaping": True,
#                 "pulse_shaping_config": {
#                     "pulse_shape": "sine_squared",
#                     "pulse_shape_rising": True,
#                     "pulse_shape_falling": True,
#                     "sample_rate_khz": 250,
#                     "rolloff_time_us": 100
#                 }
#             }
#         }
#     ]
