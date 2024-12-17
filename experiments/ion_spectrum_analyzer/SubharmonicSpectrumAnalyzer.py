import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)

from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper
import LAX_exp.experiments.eggs_heating.EGGSHeatingRDX as EGGSHeatingRDX


class SubharmonicSpectrumAnalyzer(EGGSHeatingRDX.EGGSHeatingRDX):
    """
    Experiment: Subharmonic Spectrum Analyzer

    todo: document
    """
    name = 'Subharmonic Spectrum Analyzer'
    kernel_invariants = {
        'config_eggs_heating_list', 'freq_sideband_readout_ftw_list', 'time_readout_mu_list', 'att_eggs_heating_mu',
        'freq_eggs_carrier_hz_list', 'freq_eggs_secular_hz_list',
        'phase_eggs_heating_rsb_turns_list', 'phase_eggs_heating_ch1_turns_list', 'waveform_index_to_phase_rsb_turns',
        'num_configs',
        # EGGS/phaser related
        'waveform_index_to_pulseshaper_vals',
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence',
        # subharmonic specials
        'freq_global_offset_hz', 'freq_carrier_0_offset_hz', 'freq_carrier_1_offset_hz', 'ampl_eggs_heating_carrier_pct'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",        NumberValue(default=10, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",   BooleanValue(default=True))
        self.setattr_argument("sub_repetitions",    NumberValue(default=1, precision=0, step=1, min=1, max=500))

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandcool_subsequence =     SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # EGGS RF
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([86.]),
                                                                                    CenterScan(83.20175, 0.05, 0.0005, randomize=True),
                                                                                ],
                                                                                global_min=0.005, global_max=4800, global_step=1,
                                                                                unit="MHz", scale=1, precision=6
                                                                            ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([1276.15]),
                                                                                    CenterScan(777.5, 4, 0.5, randomize=True),
                                                                                    ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, precision=3
                                                                            ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_readout_us_list",                       Scannable(
                                                                                default=[
                                                                                    ExplicitScan([120.5]),
                                                                                    RangeScan(0, 1500, 100, randomize=True),
                                                                                ],
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, precision=5
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("time_eggs_heating_us",                       NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 3, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, precision=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_ch1_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 21, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, precision=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",               NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("att_eggs_heating_db",                NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",          NumberValue(default=2., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",          NumberValue(default=2., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1500, precision=0, step=100, min=1, max=5000), group='EGGS_Heating.pulse_shaping')

        # EGGS RF - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",      BooleanValue(default=True), group='EGGS_Heating.waveform.psk')
        self.setattr_argument("num_psk_phase_shifts",           NumberValue(default=4, precision=0, step=10, min=1, max=200), group='EGGS_Heating.waveform.psk')

        # subharmonic spectrum analyzer - extras
        self.setattr_argument("freq_global_offset_mhz",                 NumberValue(default=2., precision=6, step=1., min=-10., max=10.), group=self.name)
        self.setattr_argument("freq_subharmonic_carrier_0_offset_khz",  NumberValue(default=0.1, precision=3, step=1., min=-10000., max=10000.), group=self.name)
        self.setattr_argument("freq_subharmonic_carrier_1_offset_khz",  NumberValue(default=-0.1, precision=3, step=1., min=-10000., max=10000.), group=self.name)
        self.setattr_argument("ampl_subharmonic_carrier_0_pct",         NumberValue(default=0.625, precision=2, step=10, min=0.0, max=99), group=self.name)
        self.setattr_argument("ampl_subharmonic_carrier_1_pct",         NumberValue(default=0.625, precision=2, step=10, min=0.0, max=99), group=self.name)
        self.setattr_argument("phase_subharmonic_carrier_0_turns",      NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group=self.name)
        self.setattr_argument("phase_subharmonic_carrier_1_turns",      NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group=self.name)
        self.setattr_argument("phase_oscillators_ch1_offset_turns",     PYONValue([0., 0., 0.5, 0.115, 0.]), group=self.name)

        self.setattr_argument("phase_subharmonic_carrier_0_psk_turns",  PYONValue([0., 0.5, 0., 0.5, 0.]), group=self.name)
        self.setattr_argument("phase_subharmonic_carrier_1_psk_turns",  PYONValue([0., 0.5, 0., 0.5, 0.]), group=self.name)


        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)


    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_oscillators_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, np.array(self.phase_oscillators_ch1_offset_turns))

        # convert freqs to Hz
        self.freq_global_offset_hz =    self.freq_global_offset_mhz * MHz
        self.freq_carrier_0_offset_hz = self.freq_subharmonic_carrier_0_offset_khz * kHz
        self.freq_carrier_1_offset_hz = self.freq_subharmonic_carrier_1_offset_khz * kHz

        # check that phaser oscillator frequencies are valid
        # ensure phaser amplitudes sum to less than 100%
        max_osc_freq_hz = max([
            np.max(list(self.freq_eggs_heating_secular_khz_list)) * kHz,
            self.freq_carrier_0_offset_hz,
            self.freq_carrier_1_offset_hz
        ]) + (self.freq_global_offset_mhz * MHz)
        min_osc_freq_hz = max([
            np.min(list(self.freq_eggs_heating_secular_khz_list)) * kHz,
            self.freq_carrier_0_offset_hz,
            self.freq_carrier_1_offset_hz
        ]) + (self.freq_global_offset_mhz * MHz)
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # call parent prepare method
        # note: need to set this before calling parent prepare_experiment
        self.ampl_eggs_heating_carrier_pct = 0.
        super().prepare_experiment()

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.waveform_index_to_pulseshaper_vals =   list()      # store compiled waveforms
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

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence
        _sequence_blocks = np.zeros((num_blocks, 4, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_eggs_heating_rsb_pct
        _sequence_blocks[:, 1, 0] = self.ampl_eggs_heating_bsb_pct
        _sequence_blocks[:, 2, 0] = self.ampl_subharmonic_carrier_0_pct
        _sequence_blocks[:, 3, 0] = self.ampl_subharmonic_carrier_1_pct

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        phase_bsb_update_delay_turns = (self.freq_global_offset_hz + np.mean(self.freq_eggs_secular_hz_list)) * (self.phaser_eggs.t_sample_mu * ns)
        phase_carrier_0_update_delay_turns = (self.freq_global_offset_hz + self.freq_carrier_0_offset_hz) * (2 * self.phaser_eggs.t_sample_mu * ns)
        phase_carrier_1_update_delay_turns = (self.freq_global_offset_hz + self.freq_carrier_1_offset_hz) * (3 * self.phaser_eggs.t_sample_mu * ns)

        _sequence_blocks[:, 1, 1] += phase_bsb_update_delay_turns + self.phase_eggs_heating_bsb_turns
        _sequence_blocks[:, 2, 1] += phase_carrier_0_update_delay_turns + self.phase_subharmonic_carrier_0_turns
        _sequence_blocks[:, 3, 1] += phase_carrier_1_update_delay_turns + self.phase_subharmonic_carrier_1_turns

        # set PSK phases on BOTH carriers
        if self.enable_phase_shift_keying:
            # check that phase shifting schedule is valid
            if (
                    (type(self.phase_subharmonic_carrier_0_psk_turns) is not list) or
                    (type(self.phase_subharmonic_carrier_1_psk_turns) is not list) or
                    (len(self.phase_subharmonic_carrier_0_psk_turns) != num_blocks) or
                    (len(self.phase_subharmonic_carrier_1_psk_turns) != num_blocks)
            ):
                raise ValueError("Invalid PSK schedule. Must be list with same length as num_psk_phase_shifts+1.")

            # PSK on carrier 0
            # _sequence_blocks[::2, 2, 1] +=  0.
            # _sequence_blocks[1::2, 2, 1] += self.phase_subharmonic_carrier_0_psk_turns
            _sequence_blocks[:, 2, 1] += self.phase_subharmonic_carrier_0_psk_turns

            # PSK on carrier 1
            # _sequence_blocks[::2, 3, 1] +=  0.
            # _sequence_blocks[1::2, 3, 1] += self.phase_subharmonic_carrier_1_psk_turns
            _sequence_blocks[:, 3, 1] += self.phase_subharmonic_carrier_1_psk_turns


        # record EGGS pulse waveforms
        for i in range(len(self.phase_eggs_heating_rsb_turns_list)):
            # update sequence block with rsb phase
            _sequence_blocks[:, 0, 1] += self.phase_eggs_heating_rsb_turns_list[i]

            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_blocks
            self.spinecho_wizard.calculate_pulseshape()
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_index_to_pulseshaper_vals.append(self.spinecho_wizard.get_waveform())

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
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz - self.freq_global_offset_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz - self.freq_global_offset_hz)
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
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(self.freq_global_offset_hz - sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(self.freq_global_offset_hz - sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(self.freq_global_offset_hz + sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(self.freq_global_offset_hz + sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier 0)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(self.freq_global_offset_hz + self.freq_carrier_0_offset_hz)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(self.freq_global_offset_hz + self.freq_carrier_0_offset_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 3 (carrier 1)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[3].set_frequency(self.freq_global_offset_hz + self.freq_carrier_1_offset_hz)
            self.phaser_eggs.channel[1].oscillator[3].set_frequency(self.freq_global_offset_hz + self.freq_carrier_1_offset_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
