import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXEnvironment

import matplotlib.pyplot as plt


class SpinEchoWizard(LAXEnvironment):
    """
    Helper: Spin Echo Wizard

    Design a spin-echo type pulse sequence for phaser based on user input.
    Output format is designed for Phaser Pulse Shaper object.
    """
    name = 'Spin Echo Wizard'
    # kernel_invariants = {
    #     "_max_waveforms", "t_max_phaser_update_rate_mu",
    #     "_dma_names", "_dma_handles"
    # }


    def build(self):
        # general
        self.time_pulse_us =                500

        # pulse shaping
        self.enable_pulse_shaping =         True
        self.time_pulse_shape_rolloff_us =  100
        self.freq_pulse_shape_sample_khz =  500
        self.type_pulse_shape =             'sine_squared'
        self.pulse_shape_blocks =           True

        # spin-echo
        self.enable_delay_spinecho =        True
        self.time_delay_spinecho_us =       1000
        self.sequence_blocks =              [
            [(40., -0.2), (40., 0.2), (20., 0.)],
            [(40., -0.2), (40., 0.7), (20., -0.25)],
            [(40., -0.2), (40., 0.2), (20., 0.)],
            [(40., -0.2), (40., 0.7), (20., -0.25)]
        ]

        # field geometry
        # todo: target
        # todo: phase values

        # get relevant devices
        # self.setattr_device('core')
        # self.setattr_device('phaser_eggs')

    def prepare(self):
        """
        todo: document
        """
        # structure: create time blocks
        # num blocks + delays
        # note: block has - time, osc ampls, osc phases

        # structure: calculate general amplitude windows
        # calculate a priori
        # multiply by ampls

        # structure: fill in time blocks (create lists of times and values)
        # multiply windows by osc ampls
        # set phases

        # structure: join blocks for final waveform
        # replace delays by delay
        # stitch pulse shaping with blocks


        # note: one compile function; reuse many times for sweeps
        # note: test function that displays waveforms

        # set global variables
        self._max_waveforms =               64
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators) is (conservatively) about 1.5 MSPS (i.e. 25 sample periods))
        self.t_max_phaser_update_rate_mu =  25 * self.phaser_eggs.t_sample_mu

        # convert pulse times to mu
        self.time_pulse_mu = self.core.seconds_to_mu(self.time_pulse_us * us)
        # todo: account for 40ns delay from each osc
        self.time_delay_spinecho_us = self.core.seconds_to_mu(self.time_delay_spinecho_us * us)
        # todo: ensure both pulse and spinecho times are 40ns multiples


    # todo: make general, anyone can use
    def calculate_pulseshape(self) -> TNone:
        """
        todo: document
        """
        '''calculate num samples (i.e. x-axis)'''
        # convert build variables to units of choice
        self.time_pulse_shape_rolloff_mu =          self.core.seconds_to_mu(self.time_pulse_shape_rolloff_us * us)
        self.time_pulse_shape_sample_mu =           self.core.seconds_to_mu(1. / (self.freq_pulse_shape_sample_khz * kHz))

        # ensure pulse shaping sample interval is valid (greater than min val)
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators)
        # is (conservatively) about 1.5 MSPS (i.e. 25 sample periods)
        if self.time_pulse_shape_sample_mu < self.t_max_phaser_update_rate_mu:
            raise Exception("Error: waveform sample rate too fast.")

        # ensure pulse shaping sample interval is multiple of phaser sample rate
        # note: -2 accounts for 2x t_sample_mu delay from having to set 3 oscillators
        num_multiples = round(self.time_pulse_shape_sample_mu / self.phaser_eggs.t_sample_mu) - 2
        self.time_pulse_shape_sample_mu = np.int64(num_multiples * self.self.phaser_eggs.t_sample_mu)

        # calculate number of samples
        self.num_pulse_shape_samples = round(self.time_pulse_shape_rolloff_mu / self.time_pulse_shape_sample_mu)
        # ensure rolloff time is integer number of pulse shape samples
        self.time_pulse_shape_rolloff_mu = np.int64(self.num_pulse_shape_samples * self.time_pulse_shape_rolloff_mu)

        # create x-axis value array
        self._pulse_shape_array_times_mu = np.linspace(0, self.num_pulse_shape_samples, self.num_pulse_shape_samples, dtype=np.int64)


        '''calculate pulse shape window (i.e. y-vals)'''
        if self.type_pulse_shape == 'sine_squared':
            # calculate x-axis scaling factor
            scale_factor_x = np.pi / (2. * self.time_pulse_shape_rolloff_mu)
            # calculate sine squared window
            self.ampl_window_rising = np.power(np.sin(scale_factor_x * self._pulse_shape_array_times_mu), 2)
            self.ampl_window_falling = self.ampl_window_rising[::-1]

        elif self.type_pulse_shape == 'error_function':
            raise Exception('Error: error function window not implemented')
        else:
            raise Exception('Error: idk, some window problem')


    # def create_block(self) -> TNone:
    #     """
    #     todo: document
    #     """
    #     self.sequence_blocks = np.array([
    #         [[40., -0.2], [40., 0.2], [20., 0.]],
    #         [[20., -0.2], [20., 0.7], [20., -0.25]],
    #         [[40., -0.2], [40., 0.2], [20., 0.]],
    #         [[20., -0.2], [20., 0.7], [20., -0.25]]
    #     ])
    #
    #     # todo: set up blocks
    #     # todo: don't bother doing any pulse shaping if enable_ps is OFF
    #     # todo: ind pulseshape OR delay means do fall and rise for EVERY block
    #     # todo: if 0 ampl between two blocks, then pulseshape
    #
    #     # todo: only do transitions if pulse_shape_blocks AND enable_delay are off
    #     # todo: get transitions between each block
    #     # ampl_start, ampl_stop, phase_st
    #     # self.transitions = self.sequence_blocks[:, :, 0]
    #     num_blocks, num_oscs, _ = np.shape(self.sequence_blocks)
    #     self.transitions0 = np.concatenate(([np.zeros(num_oscs)], self.sequence_blocks[:, :, 0], [np.zeros(num_oscs)]), axis=0)
    #     self.transitions10 = self.transitions0[1:, :]
    #     self.transitions11 = self.transitions0[:-1, :]
    #
    #     self.transitions2 = []
    #     for i in range(num_blocks + 1):
    #         thde = []
    #         for j in range(num_oscs):
    #             thde.append((self.transitions10[i, j], self.transitions11[i, j]))
    #
    #         self.transitions2.append(thde)
    #
    #     # todo: need ampl changes,
    #
    #
    #     # todo: think in terms of transitions - normally prev ampl to next ampl
    #     # todo: but if pulse_shape_blocks is on, then should go to 0 (i.e. still 1 transition block, but transition shape itself is different
    #
    #     # todo: actually, maybe not - way to do it is that each block has a pulse_rise, main, pulse_fall
    #
    #     # todo: should still put in a set for main blocks
    #     # todo: make delays coded as main blocks

    def compile_waveform(self):
        """
        todo: document
        record waveform as DMA sequence
        """
        self.sequence_blocks = np.array([
            [[40., -0.2], [40., 0.2], [20., 0.]],
            [[20., -0.2], [20., 0.7], [20., -0.25]],
            [[40., -0.2], [40., 0.2], [20., 0.]],
            [[20., -0.2], [20., 0.7], [20., -0.25]]
        ])

        # get shapes of blocks
        num_blocks, num_oscs, _ = np.shape(self.sequence_blocks)

        # create temporary holder objects
        _ampl_arrs =    np.zeros((num_oscs, 0), dtype=np.float64)
        _phas_arrs =    np.zeros((num_oscs, 0), dtype=np.float64)
        _time_arr =     np.array(0, dtype=np.int64)


        ###BEGIN COMPILATION - CONTIGUOUS SEQUENCES###
        if self.pulse_shape_blocks is False:

            # implement pulse shaping: rising
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


            # process blocks
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


            # implement pulse shaping: falling
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


        ###BEGIN COMPILATION - NON-CONTIGUOUS SEQUENCES###
        else:

            # process blocks
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
                        wind_phas_tmp = np.ones(len(wind_ampl_tmp)) * block_osc_vals[1]
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
                        wind_phas_tmp = np.ones(len(wind_ampl_tmp)) * block_osc_vals[1]
                        _phas_arrs[idx_osc] = np.append(_phas_arrs[idx_osc], wind_phas_tmp)

                    # do timing
                    time_tmp = np.ones(len(self.ampl_window_falling), dtype=np.int64) * self.time_pulse_shape_sample_mu
                    _time_arr = np.append(_time_arr, time_tmp)


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


    def display_waveform(self):
        """
        todo: document
        """
        pass


