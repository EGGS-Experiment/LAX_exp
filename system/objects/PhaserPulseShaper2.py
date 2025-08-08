from artiq.experiment import *
from numpy import array, ndarray, int32, int64

from LAX_exp.extensions import *
from LAX_exp.base import LAXEnvironment

PULSESHAPER_MAX_WAVEFORMS = 64


class PhaserPulseShaper2(LAXEnvironment):
    """
    Helper: Phaser Pulse Shaper 2

    Simplifies complex waveform programming on phaser via DMA and manages related DMA handles.
    Does phaser CH0 and phaser CH1 separately.
    """
    name = 'Phaser Pulse Shaper'
    kernel_invariants = {
        "t_max_phaser_update_rate_mu", "_dma_names"
    }

    def build(self, phase_offsets_turns=array([0., 0., 0., 0., 0.])):
        """
        TODO: DOCUMENT
        Arguments:
            phase_offsets_turns: todo: document
        """
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser_eggs')

        # initialize important variables
        # note: we do this here since "prepare" method is run AFTER prepare_experiment
        if (isinstance(phase_offsets_turns, list) or isinstance(phase_offsets_turns, ndarray)) and len(phase_offsets_turns) == 5:
            self.phase_offsets_turns = phase_offsets_turns
        else:
            raise ValueError("Error in PhaserPulseShaper: phase_offsets_turns must be list of length 5.")

    def prepare(self):
        """
        Prepare relevant values for waveform compilation.
        """
        # note: without touching core analyzer, max amplitude update rate for phaser (with 3 oscillators)
        # is (conservatively) about 1.5 MSPS (i.e. 25 sample periods)
        self.t_max_phaser_update_rate_mu =  25 * self.phaser_eggs.t_sample_mu

        # create data structures to allow programmatic recording & playback of DMA handles
        self._dma_names =   ['_phaser_waveform_{:d}'.format(i) for i in range(PULSESHAPER_MAX_WAVEFORMS)]
        self._dma_handles = [(0, int64(0), int32(0), False)] * PULSESHAPER_MAX_WAVEFORMS

        # store number of waveforms recorded
        self._num_waveforms = 0


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def waveform_record(self, ampl_frac_list0: TArray(TFloat, 2), ampl_frac_list1: TArray(TFloat, 2),
                        phas_turns_list0: TArray(TFloat, 2), phas_turns_list1: TArray(TFloat, 2),
                        sample_interval_mu_list: TArray(TInt64, 1)) -> TInt32:
        """
        Record waveform as DMA sequence.
        todo: finish documenting
        Arguments:
            todo
        Returns:
                (TInt32)    : the index of the recorded waveform (for later playback).
        """
        '''PREPARE INPUTS'''
        # get total lengths of arrays
        len_ampls =     len(ampl_frac_list0)
        len_phases =    len(phas_turns_list0)
        len_times =     len(sample_interval_mu_list)

        # get length of each row (i.e. number of oscs)
        num_ampl_vals = len(ampl_frac_list0[0])
        num_phas_vals = len(phas_turns_list0[0])

        # ensure all inputs have correct dimensionality
        if not ((len_ampls == len_phases) and (len_phases == len_times)):
            raise ValueError("Error: waveform arrays do not have same sizes.")

        # ensure oscillator updates have correct dimensionality
        if not (num_ampl_vals == num_phas_vals):
            raise ValueError("Error: waveform arrays do not have same sizes.")
        elif (num_ampl_vals > 5) or (num_phas_vals > 5):
            raise ValueError("Error: waveform arrays exceed number of oscillators.")

        # ensure we haven't exceeded max number of waveforms
        if self._num_waveforms > PULSESHAPER_MAX_WAVEFORMS:
            raise ValueError("Error: too many waveforms recorded.")


        '''RECORD WAVEFORMS'''
        # save waveform as dataset for posterity
        # note: do here so we don't have to break_realtime again
        self._waveform_save_dataset(self._num_waveforms,
                                    ampl_frac_list0, phas_turns_list0,
                                    ampl_frac_list1, phas_turns_list1,
                                    sample_interval_mu_list)

        # add slack for recording DMA sequences (1 ms)
        self.core.break_realtime()
        delay_mu(1000000)

        # record phaser rising pulse shape DMA sequence
        with self.core_dma.record(self._dma_names[self._num_waveforms]):
            for i in range(len_ampls):
                # set waveform values
                self._waveform_point(ampl_frac_list0[i], ampl_frac_list1[i], phas_turns_list0[i], phas_turns_list1[i])
                # set variable delay
                delay_mu(sample_interval_mu_list[i])
        # add slack after DMA recording
        self.core.break_realtime()

        # increment waveform counter and return current index
        self._num_waveforms += 1
        return self._num_waveforms - 1

    @rpc(flags={"async"})
    def _waveform_save_dataset(self, wav_idx: TInt32,
                               ampl_frac_list_0: TArray(TFloat, 2), phas_turns_list_0: TArray(TFloat, 2),
                               ampl_frac_list_1: TArray(TFloat, 2), phas_turns_list_1: TArray(TFloat, 2),
                               sample_interval_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Save waveform as dataset.
        In case a question is later raised about waveforms etc.
        """
        self.set_dataset("_phaser_wav_{:d}_CH0_ampls".format(wav_idx), ampl_frac_list_0)
        self.set_dataset("_phaser_wav_{:d}_CH0_phases".format(wav_idx), phas_turns_list_0)
        self.set_dataset("_phaser_wav_{:d}_CH1_ampls".format(wav_idx), ampl_frac_list_1)
        self.set_dataset("_phaser_wav_{:d}_CH1_phases".format(wav_idx), phas_turns_list_1)
        self.set_dataset("_phaser_wav_{:d}_times".format(wav_idx), sample_interval_mu_list)

    @kernel(flags={"fast-math"})
    def _waveform_point(self, ampl_frac_list0: TArray(TFloat, 1), ampl_frac_list1: TArray(TFloat, 1),
                        phase_turns_list0: TArray(TFloat, 1), phase_turns_list1: TArray(TFloat, 1)) -> TNone:
        """
        todo: document
        todo: note that this function can do fewer than 5 oscs
        """
        # loop over input array (guided by ampl_frac_list)
        for osc_num in range(len(ampl_frac_list0)):
            # set outputs for both phaser channels in parallel
            with parallel:
                self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list0[osc_num],
                    phase=phase_turns_list0[osc_num],
                    clr=0
                )
                self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list1[osc_num],
                    phase=phase_turns_list1[osc_num] + self.phase_offsets_turns[osc_num],
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
        # at_mu(self.phaser_eggs.get_next_frame_mu())
        self.core_dma.playback_handle(self._dma_handles[waveform_num])
