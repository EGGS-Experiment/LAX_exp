import numpy as np
from artiq.experiment import *

from LAX_exp.language import *
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

import matplotlib.pyplot as plt


class SpinEchoWizard(LAXEnvironment):
    """
    Helper: Spin Echo Wizard

    Design a spin-echo type pulse sequence for phaser based on user input.
    Output format is designed for Phaser Pulse Shaper object.
    Does not interact with the core device.
    """
    name = 'Spin Echo Wizard'
    # todo: kernel_invariants
    # kernel_invariants = {
    #     "_max_waveforms", "t_max_phaser_update_rate_mu",
    #     "_dma_names", "_dma_handles"
    # }

    def build(self):
        """
        Build SpinEchoWizard and create configuration variables.
        """
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('phaser_eggs')

        # set max hardware sample rate
        # note: max update rate should be a multiple of 5x the sample period
        # such that each oscillator is deterministically updated
        self.t_max_phaser_update_rate_mu =  5 * self.phaser_eggs.t_sample_mu
        self.num_max_phaser_samples =       190


        '''WAVEFORM CONFIGURATION ARGUMENTS'''
        # waveform block sequence
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        self.time_pulse_us = 200
        # self.sequence_blocks = np.array([
        #     [[40., -0.2], [40., 0.2], [20., 0.]],
        #     [[20., -0.2], [20., 0.7], [20., -0.25]],
        #     [[40., -0.2], [40., 0.2], [20., 0.]],
        #     [[20., -0.2], [20., 0.7], [20., -0.25]]
        # ])
        self.sequence_blocks = np.array([
            [[40., -0.0], [40., 0.0], [20., 0.]]
            # [[40., -0.0], [40., 0.0], [20., 0.]],
            # [[40., -0.0], [40., 0.0], [20., 0.]]
        ])

        # pulse shaping
        self.enable_pulse_shaping =         True
        self.pulse_shape_blocks =           False
        self.type_pulse_shape =             'sine_squared'
        self.time_pulse_shape_rolloff_us =  100
        self.freq_pulse_shape_sample_khz =  500

        # spin-echo delay
        self.enable_delay_spinecho =        False
        self.time_delay_spinecho_us =       250

    def calculate_pulseshape(self) -> TNone:
        """
        Calculate desired pulseshape based on configured options.
        """
        '''precalculate timings'''
        # convert pulse times to mu and ensure they are multiples of the phaser sample period (40ns)
        # note: 1 frame period = 4 ns/clock * 8 clock cycles * 10 words = 320ns
        self.time_pulse_mu = self.core.seconds_to_mu(self.time_pulse_us * us)
        self.time_pulse_mu -= self.time_pulse_mu % self.phaser_eggs.t_sample_mu

        self.time_delay_spinecho_mu = self.core.seconds_to_mu(self.time_delay_spinecho_us * us)
        self.time_delay_spinecho_mu -= self.time_delay_spinecho_mu % self.phaser_eggs.t_sample_mu

        '''calculate num samples (i.e. x-axis)'''
        # convert build variables to units of choice
        time_pulse_shape_rolloff_mu =  self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)
        self.time_pulse_shape_sample_mu =   self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))
        # ensure pulse shaping sample interval is valid (greater than min val)
        if self.time_pulse_shape_sample_mu < (5 * self.phaser_eggs.t_sample_mu):
            raise Exception("Error: waveform sample rate too fast.")

        # ensure pulse shaping sample interval is multiple of phaser sample rate
        # note: -2 accounts for 2x t_sample_mu delay from having to set 3 oscillators
        # note: evidently, this assumes we have exactly 3 oscillators - won't be accurate if not 3 oscs
        num_multiples = round(self.time_pulse_shape_sample_mu / self.phaser_eggs.t_sample_mu) - 2
        self.time_pulse_shape_sample_mu = np.int64(num_multiples * self.phaser_eggs.t_sample_mu)

        # calculate number of samples
        num_samples = round(time_pulse_shape_rolloff_mu / self.time_pulse_shape_sample_mu)
        # ensure rolloff time is integer number of pulse shape samples
        time_pulse_shape_rolloff_mu = np.int64(num_samples * self.time_pulse_shape_sample_mu)
        # ensure we aren't submitting too many samples to prevent overloading DMA
        if num_samples > self.num_max_phaser_samples:
            raise ValueError("Too many points in phaser waveform. Reduce sample rate or rollon time.")

        '''calculate pulse shape window (i.e. y-vals)'''
        try:
            # get corresponding pulse shape function
            pulse_shape_func = available_pulse_shapes[self.type_pulse_shape]
            self.ampl_window_rising = pulse_shape_func(
                np.arange(num_samples, dtype=np.int64) * self.time_pulse_shape_sample_mu,  # array of time-series values
                time_pulse_shape_rolloff_mu
            )

            # falling pulse shape should be simple reverse of rising window
            self.ampl_window_falling = self.ampl_window_rising[::-1]
        except KeyError:
            raise ValueError('Invalid pulse shape: {}'.format(self.type_pulse_shape))

    def compile_waveform(self) -> TNone:
        """
        Compile configured waveform for use with PhaserPulseShaper.
        todo: why store result? why not just return?
        """
        # get shapes of blocks
        num_blocks, num_oscs, _ = np.shape(self.sequence_blocks)

        # create temporary holder objects
        _ampl_arrs =    list(np.zeros((num_oscs, 0), dtype=np.float64))
        _phas_arrs =    list(np.zeros((num_oscs, 0), dtype=np.float64))
        _time_arr =     np.zeros(0, dtype=np.int64)

        # determine whether constituent blocks in the sequence are contiguous, or non-contiguous
        # contiguous: only a single pulse-shape is applied; namely to the first and last blocks of the sequence
        # non-contiguous: all constituent blocks have a pulse shape applied to them
        _contiguous_sequence = ((self.pulse_shape_blocks is False) or
                                ((self.enable_delay_spinecho is False) and (self.enable_pulse_shaping is False)))


        ###BEGIN COMPILATION - CONTIGUOUS SEQUENCES###
        if _contiguous_sequence:

            # step 1: implement pulse shape (rising)
            # note: the rising pulse shape assumes the ampl and phase of the first block in the sequence
            if self.enable_pulse_shaping is True:
                # pulse shape each oscillator individually
                for idx_osc in range(num_oscs):
                    # get block vals
                    block_osc_vals = self.sequence_blocks[0, idx_osc]

                    # do ampls
                    wind_ampl_tmp = self.ampl_window_rising * block_osc_vals[0]
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)

                    # do phases
                    wind_phas_tmp = np.ones(len(wind_ampl_tmp)) * block_osc_vals[1]
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)

                # do timing
                time_tmp = np.ones(len(self.ampl_window_rising), dtype=np.int64) * self.time_pulse_shape_sample_mu
                _time_arr = np.append(_time_arr, time_tmp)


            # step 2: process blocks
            for idx_block in range(num_blocks):
                # set values for each oscillator individually
                for idx_osc in range(num_oscs):
                    # get block vals
                    block_osc_vals = self.sequence_blocks[idx_block, idx_osc]

                    # do ampl, phase
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], block_osc_vals[0])
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], block_osc_vals[1])

                # do timing
                _time_arr = np.append(_time_arr, self.time_pulse_mu)


            # step 3: implement pulse shaping: falling
            # note: the falling pulse shape assumes the ampl and phase of the last block in the sequence
            if self.enable_pulse_shaping is True:
                # pulse shape each oscillator individually
                for idx_osc in range(num_oscs):
                    # get block vals
                    block_osc_vals = self.sequence_blocks[-1, idx_osc]

                    # do ampls
                    wind_ampl_tmp = self.ampl_window_falling * block_osc_vals[0]
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)

                    # do phases
                    wind_phas_tmp = np.ones(len(wind_ampl_tmp)) * block_osc_vals[1]
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)

                # do timing
                time_tmp = np.ones(len(self.ampl_window_falling), dtype=np.int64) * self.time_pulse_shape_sample_mu
                _time_arr = np.append(_time_arr, time_tmp)


            # step 4: ensure pulse turns off
            for idx_osc in range(num_oscs):
                _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], 0.)
                _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], 0.)
            # set minimum delay time after turnoff
            _time_arr = np.append(_time_arr, self.phaser_eggs.t_sample_mu)


        ###BEGIN COMPILATION - NON-CONTIGUOUS SEQUENCES###
        else:

            # process each block separately
            for idx_block in range(num_blocks):
                # get block vals
                block_vals = self.sequence_blocks[idx_block]

                # implement pulse shaping: rising
                if self.enable_pulse_shaping is True:
                    # pulse shape each oscillator individually
                    for idx_osc in range(num_oscs):
                        # get block vals
                        block_osc_vals = block_vals[idx_osc]

                        # do ampls
                        wind_ampl_tmp = self.ampl_window_rising * block_osc_vals[0]
                        _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)

                        # do phases
                        wind_phas_tmp = np.ones(len(self.ampl_window_rising)) * block_osc_vals[1]
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)

                    # do timing
                    time_tmp = np.ones(len(self.ampl_window_rising), dtype=np.int64) * self.time_pulse_shape_sample_mu
                    _time_arr = np.append(_time_arr, time_tmp)


                # process block - main pulse
                for idx_osc in range(num_oscs):
                    # get block vals
                    block_osc_vals = block_vals[idx_osc]

                    # do ampl, phase
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], block_osc_vals[0])
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], block_osc_vals[1])

                # do timing
                _time_arr = np.append(_time_arr, self.time_pulse_mu)


                # implement pulse shaping: falling
                if self.enable_pulse_shaping is True:
                    # pulse shape each oscillator individually
                    for idx_osc in range(num_oscs):
                        # get block vals
                        block_osc_vals = block_vals[idx_osc]

                        # do ampls
                        wind_ampl_tmp = self.ampl_window_falling * block_osc_vals[0]
                        _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], wind_ampl_tmp)

                        # do phases
                        wind_phas_tmp = np.ones(len(self.ampl_window_falling)) * block_osc_vals[1]
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)

                    # do timing
                    time_tmp = np.ones(len(self.ampl_window_falling), dtype=np.int64) * self.time_pulse_shape_sample_mu
                    _time_arr = np.append(_time_arr, time_tmp)


                # ensure pulse turns off
                for idx_osc in range(num_oscs):
                    _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], 0.)
                    _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], 0.)
                # set minimum delay time
                _time_arr = np.append(_time_arr, self.phaser_eggs.t_sample_mu)


                # implement spin-echo delay between blocks
                if self.enable_delay_spinecho and (idx_block < (num_blocks - 1)):
                    # clear ampl and phase for each osc
                    for idx_osc in range(num_oscs):
                        _ampl_arrs[idx_osc] = np.append(_ampl_arrs[idx_osc], 0.)
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], 0.)

                    # set delay time
                    _time_arr = np.append(_time_arr, self.time_delay_spinecho_mu)

        # store results
        self._ampl_tmp_arr = np.array(_ampl_arrs)
        self._phas_tmp_arr = np.array(_phas_arrs)
        self._time_tmp_arr = np.array(_time_arr)

    def get_waveform(self):
        """
        Retrieve the designed waveform.
        todo: annotate return types
        """
        # note: need to convert ampls from pct to frac
        wav_data_ampl = self._ampl_tmp_arr.transpose() / 100.
        wav_data_phas = self._phas_tmp_arr.transpose()
        wav_data_time = self._time_tmp_arr
        return wav_data_ampl, wav_data_phas, wav_data_time

    def display_waveform(self):
        """
        Display waveform phase and amplitude of each oscillator.
        todo: document
        """
        try:
            # get sequence shape
            num_blocks, num_oscs, _ = np.shape(self.sequence_blocks)

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
