import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM, ReadoutAdaptive
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper

import math
from itertools import product
from collections.abc import Iterable
from artiq.coredevice import ad9910


class SuperDuperResolutionCharacteristicReconstruction(LAXExperiment, Experiment):
    """
    Experiment: Super Duper Resolution Characteristic Reconstruction

    Characteristic Function Reconstruction (for Wigner tomography) of the Super Duper Resolution technique.
    """
    name = 'SuperRes Char Read'
    kernel_invariants = {
        # hardware objects
        'singlepass0', 'singlepass1',

        # hardware values - superresolution
        'att_eggs_heating_mu', 'freq_superresolution_osc_base_hz_list', 'freq_phaser_center_hz', 'pulseshaper_vals',
        'time_superresolution_stop_mu',

        # hardware values - characteristic readout - default
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'freq_sigmax_ftw', 'ampl_sigmax_asf', 'time_sigmax_mu', 'att_reg_sigmax',

        # hardware values - characteristic readout - characteristic actual
        'ampls_cat_asf', 'phase_characteristic_axis_pow', 'phases_pulse4_cat_pow',
        'phases_pulse4_cat_update_dir', 'freq_cat_center_ftw', 'freq_cat_secular_ftw', 'att_reg_bichromatic',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',
        # configs
        'profile_729_SBC', 'profile_729_target', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # allocate relevant beam profiles
        self.profile_729_SBC =      1
        self.profile_729_target =   6

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence =       RescueIon(self, input_type="probability")

        # set arguments for different experiment bases
        self._build_arguments_superresolution()
        self._build_arguments_characteristic_reconstruction()

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('phaser_eggs')
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.singlepass1 = self.get_device("urukul0_ch2")

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def _build_arguments_superresolution(self):
        """
        Set specific arguments for superresolution.
        """
        # superresolution - general config
        self.setattr_argument("enable_phaser", BooleanValue(default=True))
        self.setattr_argument("enable_cutoff", BooleanValue(default=False))

        # superresolution - configurable freq & sweeps
        self.setattr_argument("freq_eggs_heating_carrier_mhz", NumberValue(default=86.0, precision=6, step=0.001, min=1., max=200., unit="MHz", scale=1.),
                              group="{}.freq_phase_sweep".format("SDR"))
        self.setattr_argument("phase_eggs_heating_ch1_turns", NumberValue(default=0.2, precision=5, step=0.1, min=-1., max=1., unit="turns", scale=1.),
                              group="{}.freq_phase_sweep".format("SDR"))

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format("SDR"))

        # EGGS RF - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",  BooleanValue(default=False), group="{}.psk".format("SDR"))
        self.setattr_argument("num_psk_phase_shifts",       NumberValue(default=1, precision=0, step=10, min=1, max=200), group="{}.psk".format("SDR"))
        self.setattr_argument("phase_superresolution_rsb_psk_turns",    PYONValue([0., 0.5]), group="{}.psk".format("SDR"))
        self.setattr_argument("phase_superresolution_bsb_psk_turns",    PYONValue([0., 0.5]), group="{}.psk".format("SDR"))
        self.setattr_argument("phase_subharmonic_carrier_0_psk_turns",  PYONValue([0., 0.]), group="{}.psk".format("SDR"))
        self.setattr_argument("phase_subharmonic_carrier_1_psk_turns",  PYONValue([0., 0.]), group="{}.psk".format("SDR"))

        # superresolution - custom waveform specification
        self.setattr_argument("time_eggs_heating_us",   NumberValue(default=200, precision=2, step=500, min=0.04, max=10000000, unit="us", scale=1.),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("att_eggs_heating_db",    NumberValue(default=17., precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=2., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("freq_superresolution_osc_khz_list",  PYONValue([-702.6, 702.6, 0., 0.]),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("ampl_superresolution_osc_frac_list", PYONValue([0., 40., 12., 0.]),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("phase_superresolution_osc_turns_list",   PYONValue([0., 0., 0.5, 0.5]),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("phase_oscillators_ch1_offset_turns",     PYONValue([0., 0., 0.5, 0.5, 0.]),
                              group="{}.waveform".format("SDR"))

        # superresolution - characteristic function special
        self.setattr_argument("time_superresolution_stop_us", NumberValue(default=500, precision=2, step=500, min=0.04, max=1000000, unit="us", scale=1.),
                              group="{}.waveform".format("SDR"))

    def _build_arguments_characteristic_reconstruction(self):
        """
        Set specific arguments for Characteristic Reconstruction.
        """
        # defaults - beam values
        self.max_ampl_singlepass_pct, self.min_att_singlepass_db = (60., 7.)
        self.setattr_argument("freq_singlepass_default_mhz_list",   PYONValue([120.339, 120.339]), group='defaults.beams', tooltip="[rsb_mhz, bsb_mhz]")
        self.setattr_argument("ampl_singlepass_default_pct_list",   PYONValue([50., 0.01]), group='defaults.beams', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("att_singlepass_default_db_list",     PYONValue([7., 7.]), group='defaults.beams', tooltip="[rsb_db, bsb_db]")
        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, unit="%", scale=1.),
                              group="defaults.beams")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, unit="dB", scale=1.),
                              group="defaults.beams")

        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.0962, precision=6, step=1, min=50., max=400., unit="MHz", scale=1.),
                              group="defaults.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, unit="%", scale=1.),
                              group="defaults.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, unit="dB", scale=1.),
                              group="defaults.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=1.49, precision=3, step=5, min=0.1, max=10000, unit="us", scale=1.),
                              group="defaults.sigmax")

        # defaults - bichromatic
        self.setattr_argument("freq_cat_center_mhz",    NumberValue(default=101.0962, precision=6, step=0.001, min=50., max=400., unit="MHz", scale=1.),
                              group="defaults.bichromatic")
        self.setattr_argument("freq_cat_secular_khz",   NumberValue(default=702.6, precision=4, step=1, min=0., max=4e5, unit="kHz", scale=1.),
                              group="defaults.bichromatic")
        self.setattr_argument("ampls_cat_pct",  PYONValue([50., 54.]), group='defaults.bichromatic', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("atts_cat_db",    PYONValue([13., 13.]), group='defaults.bichromatic', tooltip="[rsb_db, bsb_db]")

        # characteristic readout: configuration
        self.setattr_argument("characteristic_axis", EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'), group='characteristic',
                              tooltip="Selects the real/imag component of the characteristic function by either applying a sigma_x operation (Imag), or not (Real)."
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_characteristic_axis_turns",  NumberValue(default=0.125, precision=5, step=0.1, min=-1.0, max=1.0, unit='turns', scale=1.),
                              group='characteristic',
                              tooltip="Sets the relative phase of the sigma_x operation used to define the real/imag axis of the characteristic function.")
        self.setattr_argument("phases_pulse4_cat_turns",    PYONValue([0., 0.]), group='characteristic', tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_pulse4_cat_phase",    EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group="characteristic",
                              tooltip="Configures how the phases of the bichromatic tones are adjusted to measure the characteristic function.")
        self.setattr_argument("characteristic_readout_sweep", EnumerationValue(['Grid', 'Phase Sweep'], default='Grid'), group='characteristic',
                              tooltip="Choose sweep type when reading out the characteristic function."
                                      "'Grid' option reads out on a 2D rectangular grid."
                                      "'Phase Sweep' option reads out on a 1D phase sweep at a single time radius.")

        # characteristic readout: readout grid
        self.setattr_argument("time_pulse4_cat_x_us_list",  Scannable(
                                                            default=[
                                                                RangeScan(-45, 45, 15, randomize=True),
                                                                ExplicitScan([50]),
                                                            ],
                                                            global_min=-100000, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5), group="characteristic.grid")
        self.setattr_argument("time_pulse4_cat_y_us_list", Scannable(
                                                            default=[
                                                                RangeScan(-45, 45, 15, randomize=True),
                                                                ExplicitScan([50]),
                                                            ],
                                                            global_min=-100000, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5), group="characteristic.grid")

        # characteristic readout: single-radius measurement (alternative to 2D grid)
        self.setattr_argument("time_char_phase_sweep_us",  NumberValue(default=10., precision=3, step=5, min=0.1, max=10000, unit="us", scale=1.),
                              group="characteristic.ph_sweep")
        self.setattr_argument("phases_char_phase_sweep_turns", Scannable(
                                                                default=[
                                                                    RangeScan(-1., 1., 42, randomize=True),
                                                                    ExplicitScan([0.]),
                                                                ],
                                                                global_min=-1., global_max=1., global_step=0.1,
                                                                unit="turns", scale=1., precision=5), group="characteristic.ph_sweep")

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE & VALIDATE INPUTS'''
        self._prepare_argument_checks()

        # prepare values for different experiment bases
        self._prepare_superresolution()
        self._prepare_characteristic_reconstruction()

        # configure waveform via pulse shaper & spin echo wizard
        # note: does this have to happen here? or can we put it inside _prepare_superres?
        self._prepare_waveform()

    def _prepare_superresolution(self):
        """
        Prepare values for superresolution.
        """
        '''SUBSEQUENCE PARAMETERS'''
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_oscillators_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, np.array(self.phase_oscillators_ch1_offset_turns))

        '''HARDWARE VALUES - CONFIG'''
        self.att_eggs_heating_mu = att_to_mu(self.att_eggs_heating_db * dB)
        self.freq_phaser_center_hz = (self.freq_eggs_heating_carrier_mhz * MHz  - self.phaser_eggs.freq_center_hz -
                                      self.freq_global_offset_mhz * MHz)
        self.freq_superresolution_osc_base_hz_list = (np.array(self.freq_superresolution_osc_khz_list) * kHz +
                                                      self.freq_global_offset_mhz * MHz)

        # superresolution/characteristic reconstruction special
        self.time_superresolution_stop_mu = self.core.seconds_to_mu(self.time_superresolution_stop_us * us)
        # ensure pre-emptive stop time is multiple of phaser frame period
        self.time_superresolution_stop_mu = ((self.time_superresolution_stop_mu // self.phaser_eggs.t_frame_mu) *
                                             self.phaser_eggs.t_frame_mu)

    def _prepare_characteristic_reconstruction(self):
        """
        Prepare values for Characteristic Reconstruction.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS - BEAMS
        '''
        # default parameters
        self.freq_singlepass_default_ftw_list = [self.singlepass0.frequency_to_ftw(freq_mhz * MHz)
                                                 for freq_mhz in self.freq_singlepass_default_mhz_list]
        self.ampl_singlepass_default_asf_list = [self.singlepass0.amplitude_to_asf(ampl_asf / 100.)
                                                 for ampl_asf in self.ampl_singlepass_default_pct_list]
        self.ampl_doublepass_default_asf =  self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)

        # sigma_x waveforms
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # bichromatic waveforms
        self.freq_cat_center_ftw =  self.qubit.frequency_to_ftw(self.freq_cat_center_mhz * MHz)
        self.freq_cat_secular_ftw = self.qubit.frequency_to_ftw(self.freq_cat_secular_khz * kHz)
        self.ampls_cat_asf =    np.array([self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=np.int32)

        # create attenuation registers
        atts_cat_mu = [att_to_mu(att_db * dB) for att_db in self.atts_cat_db]
        self.att_singlepass_default_mu_list = [att_to_mu(att_db * dB)
                                               for att_db in self.att_singlepass_default_db_list]
        self.att_reg_sigmax = 0x00000000 | (
                (att_to_mu(self.att_sigmax_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_singlepass_default_mu_list[0] << ((self.singlepass0.chip_select - 4) * 8)) |
                (self.att_singlepass_default_mu_list[1] << ((self.singlepass1.chip_select - 4) * 8))
        )
        self.att_reg_bichromatic = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (atts_cat_mu[0] << ((self.singlepass0.chip_select - 4) * 8)) |
                (atts_cat_mu[1] << ((self.singlepass1.chip_select - 4) * 8))
        )

        '''
        CONVERT VALUES TO MACHINE UNITS - PULSES/TIMINGS
        '''
        # define real/imag part of characteristic axis to read out (via a sigma_x pi/2 pulse)
        self.phase_characteristic_axis_pow =  self.qubit.turns_to_pow(self.phase_characteristic_axis_turns)
        # configure whether real/imag/both axes of the characteristic function are to be measured
        if self.characteristic_axis == "Real":          characteristic_axis_list = [False]
        elif self.characteristic_axis == "Imaginary":   characteristic_axis_list = [True]
        elif self.characteristic_axis == "Both":        characteristic_axis_list = [True, False]

        # define relative axis for characteristic readout
        self.phases_pulse4_cat_pow = np.array([self.singlepass0.turns_to_pow(phas_pow)
                                               for phas_pow in self.phases_pulse4_cat_turns], dtype=np.int32)
        if self.target_pulse4_cat_phase == 'RSB':
            self.phases_pulse4_cat_update_dir = np.array([1, 0], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'BSB':
            self.phases_pulse4_cat_update_dir = np.array([0, 1], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'RSB-BSB':
            self.phases_pulse4_cat_update_dir = np.array([1, -1], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'RSB+BSB':
            self.phases_pulse4_cat_update_dir = np.array([1, 1], dtype=np.int32)

        # create sampling grid in radial coordinates for "Grid" readout type
        if self.characteristic_readout_sweep == "Grid":
            vals_pulse4_mu_pow_list = np.array([
                [
                    self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                    self.singlepass0.turns_to_pow(np.arctan2(y_us, x_us) / (2. * np.pi))
                ]
                for x_us in self.time_pulse4_cat_x_us_list
                for y_us in self.time_pulse4_cat_y_us_list
            ], dtype=np.int64)
        # create phase-sweep array at single time value for "Phase Sweep" readout type (to save time)
        elif self.characteristic_readout_sweep == "Phase Sweep":
            vals_pulse4_mu_pow_list = np.array([
                [
                    self.core.seconds_to_mu(self.time_char_phase_sweep_us * us),
                    self.singlepass0.turns_to_pow(phase_turns)
                ]
                for phase_turns in self.phases_char_phase_sweep_turns
            ], dtype=np.int64)

        # use generator to flatten list with some tuples
        def flatten(xs):
            for x in xs:
                if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
                    yield from flatten(x)
                else:
                    yield x

        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            list(flatten(vals))
            for vals in product(
                vals_pulse4_mu_pow_list, characteristic_axis_list
            )
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''SUPERRESOLUTION - ARGUMENT CHECKS'''
        # check that input amplitude/phase arrays are valid
        if isinstance(self.ampl_superresolution_osc_frac_list, list):
            if len(self.ampl_superresolution_osc_frac_list) != 4:
                raise ValueError("Error: phaser oscillator amplitude array must have length 4.")
            elif np.sum(self.ampl_superresolution_osc_frac_list) >= 100.:
                raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude array must be a list.")

        if isinstance(self.phase_superresolution_osc_turns_list, list):
            if len(self.phase_superresolution_osc_turns_list) != 4:
                raise ValueError("Error: phaser oscillator phase array must have length 4.")
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        if not isinstance(self.freq_superresolution_osc_khz_list, list):
            raise ValueError("Error: phaser oscillator frequency array must be a list.")
        elif len(self.freq_superresolution_osc_khz_list) != 4:
            raise ValueError("Error: phaser oscillator frequency array must have length 4.")
        max_osc_freq_hz = (
                max(self.freq_superresolution_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        min_osc_freq_hz = (
                max(self.freq_superresolution_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_eggs_heating_carrier_mhz * MHz)
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_eggs_heating_carrier_mhz * MHz)
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        # check that PSK schedule is valid
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (not isinstance(psk_schedule, list)) or (len(psk_schedule) != self.num_psk_phase_shifts + 1)
            for psk_schedule in (
                self.phase_superresolution_rsb_psk_turns, self.phase_superresolution_bsb_psk_turns,
                self.phase_subharmonic_carrier_0_psk_turns, self.phase_subharmonic_carrier_1_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule. Must be a list of length num_psk_phase_shifts+1.")

        # check that cutoff time happens before end of pulse - otherwise no point
        if (self.enable_phaser and self.enable_cutoff and
                (self.time_superresolution_stop_us > self.time_eggs_heating_us)):
            raise ValueError("Invalid cutoff time. Must happen before time_superresolution_stop_us.")

        '''CHARACTERISTIC RECONSTRUCTION - ARGUMENT CHECKS'''
        # ensure single pass values are safe and valid
        if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
                for ampl_pct in self.ampl_singlepass_default_pct_list)):
            raise ValueError("Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))

        if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
                for att_db in self.att_singlepass_default_db_list)):
            raise ValueError("Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizardRDX and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.pulseshaper_vals = None        # store compiled waveforms from pulseshaper
        self.pulseshaper_id =   np.int32(0) # store waveform ID for pulseshaper

        '''PROCESS ARGUMENTS INTO WAVEFORM CONFIGS'''
        # calculate block timings
        if self.enable_phase_shift_keying:
            time_block_us = self.time_eggs_heating_us / (self.num_psk_phase_shifts + 1)

            if not self.enable_cutoff:
                num_blocks = self.num_psk_phase_shifts + 1
                block_time_list_us = [time_block_us] * num_blocks

            # implement QVSA pulse cutoffs by modifying waveform itself
            else:
                num_blocks = round(np.ceil(self.time_superresolution_stop_us / time_block_us))
                block_time_list_us = [time_block_us] * num_blocks
                if self.time_superresolution_stop_us % time_block_us != 0:
                    block_time_list_us[-1] = self.time_superresolution_stop_us % time_block_us
        else:
            num_blocks = 1
            if not self.enable_cutoff:
                block_time_list_us = [self.time_eggs_heating_us]
            # implement QVSA pulse cutoffs by modifying waveform itself
            else:
                block_time_list_us = [self.time_superresolution_stop_us]


        '''PROGRAM & COMPILE WAVEFORM'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = np.zeros((num_blocks, 4, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = np.array(self.ampl_superresolution_osc_frac_list)

        # set oscillator phases and account for oscillator update delays
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_s_list = (self.core.mu_to_seconds(self.phaser_eggs.t_sample_mu)) * np.array([0, 1, 2, 2])
        _osc_vals_blocks[:, :, 1] += (np.array(self.phase_superresolution_osc_turns_list) +
                                      self.freq_superresolution_osc_base_hz_list * t_update_delay_s_list)

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            _osc_vals_blocks[:, 0, 1] += self.phase_superresolution_rsb_psk_turns[:num_blocks]
            _osc_vals_blocks[:, 1, 1] += self.phase_superresolution_bsb_psk_turns[:num_blocks]
            _osc_vals_blocks[:, 2, 1] += self.phase_subharmonic_carrier_0_psk_turns[:num_blocks]
            _osc_vals_blocks[:, 3, 1] += self.phase_subharmonic_carrier_1_psk_turns[:num_blocks]

        # specify sequence as a dict of blocks, where each block is a dict
        _sequence_blocks = [
            {
                "oscillator_parameters": _osc_vals_blocks[i],
                "config": {
                    "time_us": block_time_list_us[i],
                    "pulse_shaping": self.enable_pulse_shaping and ((i == 0) or (i == num_blocks - 1)),
                    "pulse_shaping_config": {
                        "pulse_shape":          self.type_pulse_shape,
                        "pulse_shape_rising":   self.enable_pulse_shaping and (i == 0),
                        "pulse_shape_falling":  self.enable_pulse_shaping and (i == num_blocks - 1) and
                                                not self.enable_cutoff,
                        "sample_rate_khz":      self.freq_pulse_shape_sample_khz,
                        "rolloff_time_us":      self.time_pulse_shape_rolloff_us
                    }
                }
            } for i in range(num_blocks)
        ]

        # create QVSA waveform and store data in a holder
        self.pulseshaper_vals = self.spinecho_wizard.compile_waveform(_sequence_blocks)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        """
        Initialize experiment hardware immediately before kernel.
        """
        '''CHARACTERISTIC FUNCTION INITIALIZATION'''
        # ensure phase_autoclear disabled on all beams to prevent phase accumulator reset
        # enable RAM mode and clear DDS phase accumulator
        self.qubit.set_cfr1()
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # set up singlepass AOMs to default values (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass_default_ftw_list[0],
                                      asf=self.ampl_singlepass_default_asf_list[0],
                                      profile=i)
            self.singlepass1.set_mu(self.freq_singlepass_default_ftw_list[1],
                                      asf=self.ampl_singlepass_default_asf_list[1],
                                      profile=i)
            delay_mu(8000) # 8us

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000) # 25us

        '''SUPERRESOLUTION INITIALIZATION'''
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.core.break_realtime()

        ### PHASER INITIALIZATION ###
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # load waveform DMA handles
        self.pulse_shaper.waveform_load()
        delay_mu(250000) # 250us

        # create local variables for use
        _loop_iter = 0  # used to check_termination more frequently
        ion_state = (-1, 0, np.int64(0))    # store ion state


        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                time_char_read_mu =         config_vals[0]
                phase_char_read_pow =       np.int32(config_vals[1])
                characteristic_axis_bool =  bool(config_vals[2])

                # set phases for characteristic readout
                phase_char_read_pow_list = [
                    self.phases_pulse4_cat_pow[0] + self.phases_pulse4_cat_update_dir[0] * phase_char_read_pow,
                    self.phases_pulse4_cat_pow[1] + self.phases_pulse4_cat_update_dir[1] * phase_char_read_pow,
                ]
                self.core.break_realtime()

                # set phaser frequency
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC) - freq_eggs_heating_carrier_hz - phaser_eggs.freq_center_hz - freq_global_offset_hz
                    self.freq_phaser_center_hz,
                    # oscillator frequencies
                    [self.freq_superresolution_osc_base_hz_list[0], self.freq_superresolution_osc_base_hz_list[1],
                     self.freq_superresolution_osc_base_hz_list[2], self.freq_superresolution_osc_base_hz_list[3], 0.],
                    self.phase_eggs_heating_ch1_turns
                )
                delay_mu(20000)

                '''MOTIONAL STATE PREPARATION'''
                # get current time
                t_phaser_start_mu = now_mu()

                # initialize ion in S-1/2 state & sideband cool to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # set target profile to ensure we run correctly (for characteristic reconstruction)
                self.qubit.set_profile(self.profile_729_target)
                self.qubit.cpld.io_update.pulse_mu(8)

                # apply superresolution interaction
                if self.enable_phaser:
                    t_phaser_start_mu = self.phaser_run(self.pulseshaper_id)

                '''CHARACTERISTIC RECONSTRUCTION'''
                # sigma_x to select axis (does dummy if characteristic_axis_bool is False)
                self.pulse_sigmax(t_phaser_start_mu, self.phase_characteristic_axis_pow, characteristic_axis_bool)

                # bichromatic pulse + state detection
                self.pulse_bichromatic(t_phaser_start_mu, time_char_read_mu, phase_char_read_pow_list,
                                       self.freq_cat_center_ftw, self.freq_cat_secular_ftw)
                ion_state = self.readout_subsequence.run()

                '''LOOP CLEANUP'''
                # retrieve results & update dataset
                self.update_results(characteristic_axis_bool, ion_state[0],
                                    time_char_read_mu, phase_char_read_pow)
                self.core.break_realtime()

                # resuscitate ion & detect deaths
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(ion_state[0])
                self.core.break_realtime()

                # check termination more frequently in case reps are low
                if _loop_iter % 50 == 0:
                    self.check_termination()
                    self.core.break_realtime()
                _loop_iter += 1

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        """
        Clean up the experiment.
        """
        # set up singlepass AOMs to default values (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass_default_ftw_list[0],
                                    asf=self.ampl_singlepass_default_asf_list[0],
                                    profile=i)
            self.singlepass1.set_mu(self.freq_singlepass_default_ftw_list[1],
                                    asf=self.ampl_singlepass_default_asf_list[1],
                                    profile=i)
            delay_mu(8000)

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)


    '''
    HELPER FUNCTIONS - PHASER
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TInt64:
        """
        Run the main QVSA pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        Returns:
            the start time of the phaser oscillator waveform (in machine units, 64b int).
            Useful to synchronize device operation.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_eggs_heating_mu, self.att_eggs_heating_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        t_start_mu = self.phaser_eggs.get_next_frame_mu()
        at_mu(t_start_mu)
        self.pulse_shaper.waveform_playback(waveform_id)

        # EGGS - STOP
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

        # return phaser osc start time (in case others want to sync)
        return t_start_mu

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # get waveform for given sweep phase
        _wav_data_ampl, _wav_data_phas, _wav_data_time = self.pulseshaper_vals

        # record phaser pulse sequence and save returned waveform ID
        delay_mu(1000000)  # 1ms
        self.pulseshaper_id = self.pulse_shaper.waveform_record(_wav_data_ampl, _wav_data_phas, _wav_data_time)
        self.core.break_realtime()


    '''
    HELPER FUNCTIONS - CHARACTERISTIC READOUT
    '''
    @kernel(flags={"fast-math"})
    def pulse_sigmax(self, time_start_mu: TInt64 = -1, phas_pow: TInt32 = 0x0, is_real: TBool = False) -> TNone:
        """
        Run a phase-coherent sigma_x pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            phas_pow: relative phase offset for the beam.
            is_real: whether to actually run the pulse (True) or a dummy pulse (False).
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            self.freq_sigmax_ftw, asf=self.ampl_sigmax_asf, pow_=phas_pow,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0], asf=self.ampl_singlepass_default_asf_list[0], pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1], asf=self.ampl_singlepass_default_asf_list[1], pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_sigmax)

        # run pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        if is_real:
            self.qubit.on()
        else:
            self.qubit.off()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

    @kernel(flags={"fast-math"})
    def pulse_bichromatic(self, time_start_mu: TInt64, time_pulse_mu: TInt64, phas_pow_list: TList(TInt32),
                          freq_carrier_ftw: TInt32, freq_secular_ftw: TInt32) -> TNone:
        """
        Run a phase-coherent bichromatic pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            time_pulse_mu: length of pulse (in machine units).
            phas_pow_list: relative phase offset for the beams (RSB, BSB) (in pow).
            freq_carrier_ftw: carrier frequency (set by the double pass) in FTW.
            freq_secular_ftw: bichromatic separation frequency (from central frequency) in FTW.
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            freq_carrier_ftw, asf=self.ampl_doublepass_default_asf,
            pow_=0, profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0]-freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1]+freq_secular_ftw, asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_bichromatic)

        # run bichromatic pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

