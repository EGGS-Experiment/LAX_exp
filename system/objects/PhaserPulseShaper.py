import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXEnvironment


class PhaserPulseShaper(LAXEnvironment):
    """
    Helper: Phaser Pulse Shaper

    todo: document
    """
    name = 'Phaser Pulse Shaper'
    kernel_invariants = {
        "_max_waveforms", "t_max_phaser_update_rate_mu",
        "_phase_offsets_turns",
        "_dma_names"
    }


    def build(self):
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser_eggs')

    def prepare(self):
        """
        todo: document
        """
        # set global variables
        self._max_waveforms =               64
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators)
        # is (conservatively) about 1.5 MSPS (i.e. 25 sample periods)
        self.t_max_phaser_update_rate_mu =  25 * self.phaser_eggs.t_sample_mu

        # store global CH1 offsets
        self._phase_offsets_turns =         np.array([0., 0., 0., 0., 0.])

        # create data structures to allow programmatic recording & playback of DMA handles
        self._dma_names =       ['_phaser_waveform_{:d}'.format(i) for i in range(self._max_waveforms)]
        self._dma_handles =     [(0, np.int64(0), np.int32(0))] * self._max_waveforms

        # store number of waveforms recorded
        self._num_waveforms =   0


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def waveform_record(self, ampl_frac_list: TArray(TFloat, 2), phas_turns_list: TArray(TFloat, 2),
                        sample_interval_mu_list: TArray(TInt64, 1)) -> TInt32:
        """
        Record waveform as DMA sequence.
        todo: finish documenting
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
    def _waveform_point(self, ampl_frac_list: TArray(TFloat, 1), phase_turns_list: TArray(TFloat, 1)) -> TNone:
        """
        todo: document
        todo: note that this function can do fewer than 5 oscs
        """
        # loop over input array (guided by ampl_frac_list)
        for osc_num in range(len(ampl_frac_list)):
            # set outputs for both phaser channels in parallel
            # todo: account for field geometry & offsets - use phase offset addition
            with parallel:
                self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(amplitude=ampl_frac_list[osc_num],
                                                                                    phase=phase_turns_list[osc_num],
                                                                                    clr=0)
                self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(amplitude=ampl_frac_list[osc_num],
                                                                                    phase=phase_turns_list[osc_num] + self._phase_offsets_turns[osc_num],
                                                                                    clr=0)
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
        """
        # note: don't synchronize to frame - let user do this since user may also need to reset DUC
        self.core_dma.playback_handle(self._dma_handles[waveform_num])
