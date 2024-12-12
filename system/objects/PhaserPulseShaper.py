import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXEnvironment


class PhaserPulseShaper(LAXEnvironment):
    """
    Helper: Phaser Pulse Shaper

    Simplifies complex waveform programming on phaser via DMA and manages related DMA handles.
    """
    name = 'Phaser Pulse Shaper'
    kernel_invariants = {
        "_max_waveforms", "t_max_phaser_update_rate_mu",
        "_phase_offsets_turns",
        "_dma_names"
    }

    def build(self, phase_offsets_turns=np.array([0., 0., 0., 0., 0.])):
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser_eggs')

        # initialize important variables
        # note: we do this here since "prepare" method is run AFTER prepare_experiment
        if len(phase_offsets_turns) == 5:
            self._phase_offsets_turns = phase_offsets_turns
        else:
            raise Exception("Error in PhaserPulseShaper - phase_offsets_turns must have length 5: {}".format(phase_offsets_turns))

    def prepare(self):
        """
        Prepare relevant values for waveform compilation.
        """
        # set global variables
        self._max_waveforms = 64
        # note: max update rate should be a multiple of 5x the sample period
        # such that each oscillator is deterministically updated
        self.t_max_phaser_update_rate_mu = 5 * self.phaser_eggs.t_sample_mu

        # create data structures to allow programmatic recording & playback of DMA handles
        self._dma_names =   ['_phaser_waveform_{:d}'.format(i) for i in range(self._max_waveforms)]
        self._dma_handles = [(0, np.int64(0), np.int32(0), False)] * self._max_waveforms

        # store number of waveforms recorded
        self._num_waveforms = 0


    '''
    MAIN FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def waveform_record(self, ampl_frac_list: TArray(TFloat, 2), phas_turns_list: TArray(TFloat, 2),
                        sample_interval_mu_list: TArray(TInt64, 1)) -> TInt32:
        """
        Record waveform as DMA sequence.
        ampl_frac_list, phas_turns_list, and sample_interval_mu_list must have the same overall length (i.e. axis 0 length).
        ampl_frac_list and phas_turns_list must also be specified for the same number of oscillators (i.e. axis 1 length).
        Can specify updates for anywhere between [1, 5] oscillators, but the oscillator numbers must be contiguous.
        Arguments:
            ampl_frac_list: 2D array of amplitudes (fractional) for each oscillator.
                Axis 0 is list of updates for each timestamp, axis 1 is ampl for each osc.
            phas_turns_list: 2D array of phases (in turns) for each oscillator.
                Axis 0 is list of updates for each timestamp, axis 1 is phase for each osc.
            sample_interval_mu_list: 1D array of timestamps (in machine units) for each ampl/phase update.
        Returns:
                the index of the recorded waveform (for later playback).
        """
        '''PREPARE INPUTS'''
        # get total lengths of arrays
        len_ampls =     len(ampl_frac_list)
        len_phases =    len(phas_turns_list)
        len_times =     len(sample_interval_mu_list)

        # get length of each row (i.e. number of oscs)
        num_ampl_vals = len(ampl_frac_list[0])
        num_phas_vals = len(phas_turns_list[0])

        # ensure all inputs have correct dimensionality
        if not ((len_ampls == len_phases) and (len_phases == len_times)):
            raise Exception("Error: waveform arrays do not have same sizes.")

        # ensure oscillator updates have correct dimensionality
        if not (num_ampl_vals == num_phas_vals):
            raise Exception("Error: waveform arrays do not have same sizes.")
        elif (num_ampl_vals > 5) or (num_phas_vals > 5):
            raise Exception("Error: waveform arrays exceed number of oscillators.")

        # ensure we haven't exceeded max number of waveforms
        if self._num_waveforms > self._max_waveforms:
            raise Exception("Error: too many waveforms recorded.")


        '''RECORD WAVEFORMS'''
        # add slack for recording DMA sequences (1 ms)
        self.core.break_realtime()
        delay_mu(1000000)

        # record phaser rising pulse shape DMA sequence
        with self.core_dma.record(self._dma_names[self._num_waveforms]):
            for i in range(len_ampls):
                # set waveform values
                self._waveform_point(ampl_frac_list[i], phas_turns_list[i])
                # set variable delay
                delay_mu(sample_interval_mu_list[i])
        # add slack after DMA recording
        self.core.break_realtime()

        # increment waveform counter and return current index
        self._num_waveforms += 1
        return self._num_waveforms - 1

    @kernel(flags={"fast-math"})
    def _waveform_point(self, ampl_frac_list: TArray(TFloat, 1), phas_turns_list: TArray(TFloat, 1)) -> TNone:
        """
        Records update for all oscillators for a single "update" event.
        Arguments:
            ampl_frac_list: 1D array of amplitudes (fractional) for each oscillator.
            phas_turns_list: 1D array of phases (in turns) for each oscillator.
        """
        # loop over input array (guided by ampl_frac_list)
        for osc_num in range(len(ampl_frac_list)):
            # set outputs for both phaser channels in parallel
            with parallel:
                self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list[osc_num],
                    phase=phas_turns_list[osc_num],
                    clr=0
                )
                self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list[osc_num],
                    phase=phas_turns_list[osc_num] + self._phase_offsets_turns[osc_num],
                    clr=0
                )
                delay_mu(self.phaser_eggs.t_sample_mu)

    @kernel(flags={"fast-math"})
    def waveform_load(self) -> TNone:
        """
        Retrieve all waveforms from DMA.
        Must be called after all DMA sequences are recorded.
        """
        # get waveform DMA sequence handles
        for i in range(self._num_waveforms):
            self._dma_handles[i] = self.core_dma.get_handle(self._dma_names[i])
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def waveform_playback(self, waveform_num: TInt32) -> TNone:
        """
        Play back a previously recorded waveform.
        Arguments:
            waveform_num    (TInt32): The waveform index to play back.
        """
        # note: don't synchronize to frame - let user do this since user may also need to reset DUC
        self.core_dma.playback_handle(self._dma_handles[waveform_num])
