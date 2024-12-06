import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper2 import PhaserPulseShaper2

import LAX_exp.experiments.eggs_heating.EGGSHeatingRDX as EGGSHeatingRDX

# todo: create separate CH0 and CH1 attenuations
# todo: need to configure channels differently
# todo: need to use pulse shaper differently
# todo: ensure PSKing is done correctly


class EGGSHeatingRDXBiasTee(EGGSHeatingRDX.EGGSHeatingRDX):
    """
    Experiment: EGGS Heating RDX Bias Tee

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and measure ion temperature
    via sideband thermometry.
    """
    name = 'EGGS Heating RDX Bias Tee'
    kernel_invariants = {
        'config_eggs_heating_list', 'freq_sideband_readout_ftw_list', 'time_readout_mu_list', 'att_eggs_heating_mu',
        'freq_eggs_carrier_hz_list', 'freq_eggs_secular_hz_list',
        'phase_eggs_heating_rsb_turns_list', 'phase_eggs_heating_ch1_turns_list', 'waveform_index_to_phase_rsb_turns',
        'num_configs',
        # EGGS/phaser related
        'waveform_index_to_pulseshaper_vals0', 'waveform_index_to_pulseshaper_vals1',
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence'
    }

    def build_experiment(self):
        super().build_experiment()

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)
        # set correct phase delays for field geometries (0 for all b/c we're addressing separate electrode pairs)
        self.pulse_shaper = PhaserPulseShaper2(self, np.array([0., 0., 0., 0., 0.]))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.waveform_index_to_pulseshaper_vals0 =  list()      # store compiled waveforms - CH0
        self.waveform_index_to_pulseshaper_vals1 =  list()      # store compiled waveforms - CH1
        self.waveform_index_to_pulseshaper_id =     np.zeros(len(self.phase_eggs_heating_rsb_turns_list), dtype=np.int32)   # store pulseshaper waveform ID

        # set up blocks for pulse sequence
        num_blocks = 1
        if self.enable_phase_shift_keying:  num_blocks = self.num_psk_phase_shifts + 1

        # set up the spin echo wizard generally
        # note: time_pulse_us divided by num_blocks to split it equally
        self.spinecho_wizard.time_pulse_us =                self.time_eggs_heating_us / num_blocks
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           False
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        False
        self.spinecho_wizard.time_delay_spinecho_us =       250

        '''DESIGN WAVEFORM SEQUENCE - CH0'''
        # create bare waveform block sequence
        _sequence_blocks = np.zeros((num_blocks, 2, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_eggs_heating_rsb_pct
        _sequence_blocks[:, 1, 0] = self.ampl_eggs_heating_bsb_pct

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        phase_bsb_update_delay_turns = np.mean(self.freq_eggs_secular_hz_list) * (self.phaser_eggs.t_sample_mu * ns)
        _sequence_blocks[:, 1, 1] = self.phase_eggs_heating_bsb_turns + phase_bsb_update_delay_turns

        # record EGGS pulse waveforms
        for i in range(len(self.phase_eggs_heating_rsb_turns_list)):
            # update sequence block with rsb phase
            phase_rsb_turns = self.phase_eggs_heating_rsb_turns_list[i]
            _sequence_blocks[:, 0, 1] = phase_rsb_turns

            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_blocks
            self.spinecho_wizard.calculate_pulseshape()
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_index_to_pulseshaper_vals0.append(self.spinecho_wizard.get_waveform())

        '''DESIGN WAVEFORM SEQUENCE - CH1'''
        # create bare waveform block sequence
        _sequence_blocks = np.zeros((num_blocks, 2, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_eggs_heating_carrier_pct

        # set PSK phases on the carrier
        if self.enable_phase_shift_keying:
            _sequence_blocks[::2, 0, 1] =   0.
            _sequence_blocks[1::2, 0, 1] =  0.5

        # record EGGS pulse waveforms
        for i in range(len(self.phase_eggs_heating_rsb_turns_list)):
            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_blocks
            self.spinecho_wizard.calculate_pulseshape()
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_index_to_pulseshaper_vals1.append(self.spinecho_wizard.get_waveform())


    '''
    HELPER FUNCTIONS - PHASER
    '''

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each RSB phase
        for i in range(len(self.phase_eggs_heating_rsb_turns_list)):

            # get waveform for given RSB phase
            _wav_data_ampl0, _wav_data_phas0, _wav_data_time = self.waveform_index_to_pulseshaper_vals0[i]
            self.core.break_realtime()
            _wav_data_ampl1, _wav_data_phas1, _wav_data_time = self.waveform_index_to_pulseshaper_vals1[i]
            self.core.break_realtime()

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            _wav_idx = self.pulse_shaper.waveform_record(_wav_data_ampl0, _wav_data_ampl1,
                                                         _wav_data_phas0, _wav_data_phas1,
                                                         _wav_data_time)
            self.waveform_index_to_pulseshaper_id[i] = _wav_idx
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, phase_ch1_offset_turns: TFloat) -> TNone:
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the carrier frequency (in Hz).
            sideband_freq_hz        (float)     : the sideband frequency (in Hz).
            phase_ch1_offset_turns  (float)     : the phase offset for CH1 (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        phase_ch1_turns = phase_ch1_offset_turns + (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns)

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(phase_ch1_turns)
        self.phaser_eggs.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 (RSB/CH0, Carrier/CH1)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB/CH0)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
