from artiq.experiment import *
from numpy import array, int32, int64, ndarray

from LAX_exp.extensions import *
from LAX_exp.base import LAXEnvironment

# max waveforms recordable onto DMA with PhaserPulseShaper
PULSESHAPER_MAX_WAVEFORMS = 201

# relative delays between oscillators
PHASER_OSC_DELAY_NS = [0, 40e-9, 80e-9, 80e-9, 120e-9]

# indices of components inside an osc_val_block
_IDX_OSC_AMPL = 0
_IDX_OSC_PHAS = 1


class PhaserPulseShaper(LAXEnvironment):
    """
    Helper: Phaser Pulse Shaper

    Simplifies complex waveform programming on phaser via DMA and manages related DMA handles.
    """
    name = 'Phaser Pulse Shaper'
    kernel_invariants = {
        "t_max_phaser_update_rate_mu", "_phase_offsets_turns",
        "_dma_names",

        "phas_ch1_osc_turns_arr", "ampl_ch1_osc_scale_arr",
    }

    def build(self, phase_offsets_turns=array([0., 0., 0., 0., 0.])):
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser_eggs')

        # initialize important variables
        # note: we do this here since "prepare" method is run AFTER prepare_experiment
        if (isinstance(phase_offsets_turns, list) or isinstance(phase_offsets_turns, ndarray)) and len(phase_offsets_turns) == 5:
            self._phase_offsets_turns = phase_offsets_turns
        else:
            raise ValueError("Error in PhaserPulseShaper: phase_offsets_turns must be list of length 5.")

        # tmp remove - create initial copy of key values so we have them by build
        # note: max update rate should be a multiple of 5x the sample period
        # such that each oscillator is deterministically updated
        self.t_max_phaser_update_rate_mu = 5 * self.phaser_eggs.t_sample_mu

        # create data structures to allow programmatic recording & playback of DMA handles
        self._dma_names =   ['_phaser_waveform_{:d}'.format(i) for i in range(PULSESHAPER_MAX_WAVEFORMS)]
        self._dma_handles = [(0, int64(0), int32(0), False)] * PULSESHAPER_MAX_WAVEFORMS

        # store number of waveforms recorded
        self._num_waveforms = 0
        # tmp remove - create initial copy of key values so we have them by build

    def prepare(self):
        """
        Prepare relevant values for waveform compilation.
        """
        # note: max update rate should be a multiple of 5x the sample period
        # such that each oscillator is deterministically updated
        self.t_max_phaser_update_rate_mu = 5 * self.phaser_eggs.t_sample_mu

        # create data structures to allow programmatic recording & playback of DMA handles
        self._dma_names =   ['_phaser_waveform_{:d}'.format(i) for i in range(PULSESHAPER_MAX_WAVEFORMS)]
        self._dma_handles = [(0, int64(0), int32(0), False)] * PULSESHAPER_MAX_WAVEFORMS

        # store number of waveforms recorded
        self._num_waveforms = 0

        # get CH1 adjustment values per oscillator
        self.phas_ch1_osc_turns_arr = self.get_parameter('phas_ch1_osc_turns_arr', group='devices.phaser.ch1', override=False)
        if not all((
            isinstance(self.phas_ch1_osc_turns_arr, list),
            len(self.phas_ch1_osc_turns_arr) == 5,
            (isinstance(val, (int, float)) for val in self.phas_ch1_osc_turns_arr)
        )):
            raise ValueError("Invalid phas_ch1_osc_turns_arr specified in dataset manager ({:}).".format(self.phas_ch1_osc_turns_arr))

        self.ampl_ch1_osc_scale_arr = array(self.get_parameter('ampl_ch1_osc_scale_arr', group='devices.phaser.ch1', override=False))
        if not all((
            isinstance(self.ampl_ch1_osc_scale_arr, (list, ndarray)),
            len(self.ampl_ch1_osc_scale_arr) == 5,
            (isinstance(val, (int, float)) for val in self.ampl_ch1_osc_scale_arr)
        )):
            raise ValueError("Invalid ampl_ch1_osc_scale_arr specified in dataset manager ({:}).".format(self.ampl_ch1_osc_scale_arr))


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
        :param ampl_frac_list: 2D array of amplitudes (fractional) for each oscillator.
            Axis 0 is list of updates for each timestamp, axis 1 is ampl for each osc.
        :param phas_turns_list: 2D array of phases (in turns) for each oscillator.
            Axis 0 is list of updates for each timestamp, axis 1 is phase for each osc.
        :param sample_interval_mu_list: 1D array of timestamps (in machine units) for each ampl/phase update.
        :return: the index of the recorded waveform (for later playback).
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
        self._waveform_save_dataset(self._num_waveforms, ampl_frac_list, phas_turns_list, sample_interval_mu_list)

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
        self.core.break_realtime() # add slack after DMA recording

        # increment waveform counter and return waveform index
        self._num_waveforms += 1
        return self._num_waveforms - 1

    @rpc(flags={"async"})
    def _waveform_save_dataset(self, wav_idx: TInt32, ampl_frac_list: TArray(TFloat, 2),
                               phas_turns_list: TArray(TFloat, 2),
                               sample_interval_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Save waveform as dataset.
        In case a question is later raised about waveforms etc.
        :param wav_idx: todo: document
        :param ampl_frac_list: todo: document
        :param phas_turns_list: todo: document
        :param sample_interval_mu_list: todo: document
        """
        self.set_dataset("_phaser_wav_{:d}_ampls".format(wav_idx), ampl_frac_list)
        self.set_dataset("_phaser_wav_{:d}_phases".format(wav_idx), phas_turns_list)
        self.set_dataset("_phaser_wav_{:d}_times".format(wav_idx), sample_interval_mu_list)

    @kernel(flags={"fast-math"})
    def _waveform_point(self, ampl_frac_list: TArray(TFloat, 1), phas_turns_list: TArray(TFloat, 1)) -> TNone:
        """
        Records update for all oscillators for a single "update" event.
        :param ampl_frac_list: 1D array of amplitudes (fractional) for each oscillator.
        :param phas_turns_list: 1D array of phases (in turns) for each oscillator.
        """
        # loop over input array (guided by ampl_frac_list)
        for osc_num in range(len(ampl_frac_list)):
            # set outputs for both phaser channels simultaneously
            with parallel:
                self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list[osc_num],
                    phase=phas_turns_list[osc_num],
                    clr=0
                )
                self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(
                    amplitude=ampl_frac_list[osc_num] * self.ampl_ch1_osc_scale_arr[osc_num], # adjust CH1 ampls
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
        :param waveform_num: The waveform index to play back.
        """
        # note: don't synchronize to frame - let user do this since user may also need to reset DUC
        self.core_dma.playback_handle(self._dma_handles[waveform_num])
