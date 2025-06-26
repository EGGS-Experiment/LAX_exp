import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM
)

from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper

import math
from itertools import product
from collections.abc import Iterable
from artiq.coredevice import ad9910


class GKPCharacteristicReconstruction(LAXExperiment, Experiment):
    """
    Experiment: Super Duper Resolution Characteristic Reconstruction

    Characteristic Function Reconstruction (for Wigner tomography) of the Super Duper Resolution technique.
    """
    name = 'GKP Char Read'
    kernel_invariants = {
        # hardware objects
        'singlepass0', 'singlepass1',

        # hardware values - superresolution
        'att_eggs_heating_mu', 'freq_eggs_heating_carrier_hz', 'freq_superresolution_osc_base_hz_list',
        'freq_global_offset_hz', 'pulseshaper_vals', 'time_superresolution_stop_mu',
        # hardware values - characteristic readout - default
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'att_doublepass_default_mu',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'att_sigmax_mu', 'time_sigmax_mu',
        # hardware values - characteristic readout - characteristic actual
        'ampls_cat_asf', 'atts_cat_mu', 'phase_characteristic_axis_pow', 'phases_pulse4_cat_pow',
        'phases_pulse4_cat_update_dir', 'freq_cat_center_ftw', 'freq_cat_secular_ftw',

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
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # set arguments for different experiment bases
        self._build_arguments_gkp()
        self._build_arguments_characteristic_reconstruction()

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('phaser_eggs')
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.singlepass1 = self.get_device("urukul0_ch2")

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)

    def _build_arguments_gkp(self):
        """
        Set specific arguments for superresolution.
        """
        # superresolution - geeneral config
        self.setattr_argument("enable_phaser", BooleanValue(default=True))
        self.setattr_argument("enable_cutoff", BooleanValue(default=False))

        # gkp - configurable freq & sweeps
        self.setattr_argument("freq_carrier_mhz", NumberValue(default=86.0, precision=6, step=0.001, min=1., max=200.),
                              group="{}.freq_phase_sweep".format(self.name))
        self.setattr_argument("phase_ch1_turns", NumberValue(default=0.2, precision=5, step=0.1, min=-1., max=1.),
                              group="{}.freq_phase_sweep".format(self.name))

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=False), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1500, precision=0, step=100, min=1, max=5000), group='EGGS_Heating.pulse_shaping')

        # superresolution - custom waveform specification
        self.setattr_argument("time_squeezing", NumberValue(default=200, precision=2, step=500, min=0.04, max=100000000),
                              group="{}.squeezing".format(self.name))
        self.setattr_argument("freq_squeezing_khz",  NumberValue(default=702, precision=2, step=0.01, min=600, max=2000),
                              group="{}.squeezing".format(self.name))
        self.setattr_argument("ampl_squeezing_frac", NumberValue(default=10, precision=2, step=0.01, min=0.04, max=.5),
                              group="{}.squeezing".format(self.name))
        self.setattr_argument("phase_squeezing_turns", NumberValue(default=0., precision=2, step=0.01, min=0.0, max=np.pi),
                              group="{}.waveform".format(self.name))


        # self.setattr_argument("time_displacement",   NumberValue(default=200, precision=2, step=500, min=0.04, max=100000000),
        #                       group="{}.waveform".format(self.name))
        # self.setattr_argument("att_displacement_db",    NumberValue(default=17., precision=1, step=0.5, min=0, max=31.5), group="{}.waveform".format(self.name))
        # self.setattr_argument("freq_superresolution_osc_khz_list",      PYONValue([-702.2, 702.2, 0., 0.]), group="{}.waveform".format(self.name))
        # self.setattr_argument("ampl_superresolution_osc_frac_list",     PYONValue([40., 40., 10., 0.]), group="{}.waveform".format(self.name))
        # self.setattr_argument("phase_superresolution_osc_turns_list",   PYONValue([0., 0., 0., 0.]), group="{}.waveform".format(self.name))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=2., precision=6, step=1., min=-10., max=10.), group="{}.waveform".format(self.name))
        self.setattr_argument("phase_oscillators_ch1_offset_turns",     PYONValue([0., 0., 0.5, 0.5, 0.]), group="{}.waveform".format(self.name))

        # superresolution - characteristic function special
        self.setattr_argument("time_gkp_stop_us", NumberValue(default=500, precision=2, step=500, min=0.04, max=100000000),
                              group="{}.waveform".format(self.name))

    def _build_arguments_characteristic_reconstruction(self):
        """
        Set specific arguments for Characteristic Reconstruction.
        """
        # defaults - beam values
        self.max_ampl_singlepass_pct, self.min_att_singlepass_db = (58., 6.)
        self.setattr_argument("freq_singlepass_default_mhz_list",   PYONValue([120.339, 120.339]), group='defaults.beams', tooltip="[rsb_mhz, bsb_mhz]")
        self.setattr_argument("ampl_singlepass_default_pct_list",   PYONValue([50., 50.]), group='defaults.beams', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("att_singlepass_default_db_list",     PYONValue([7., 7.]), group='defaults.beams', tooltip="[rsb_db, bsb_db]")
        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.beams")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.beams")

        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.1202, precision=6, step=1, min=50., max=400.), group="defaults.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=1.55, precision=2, step=5, min=0.1, max=10000), group="defaults.sigmax")

        # defaults - bichromatic
        self.setattr_argument("freq_cat_center_mhz", NumberValue(default=101.1202, precision=6, step=0.001, min=50., max=400.), group="defaults.bichromatic")
        self.setattr_argument("freq_cat_secular_khz", NumberValue(default=702.2, precision=4, step=1, min=0., max=4e5), group="defaults.bichromatic")
        self.setattr_argument("ampls_cat_pct",  PYONValue([50., 50.]), group='defaults.bichromatic', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("atts_cat_db",    PYONValue([13., 13.]), group='defaults.bichromatic', tooltip="[rsb_db, bsb_db]")

        # characteristic readout: configure readout grid
        self.setattr_argument("characteristic_axis", EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'), group='characteristic_readout',
                              tooltip="Selects the real/imag component of the characteristic function by either applying a sigma_x operation (Imag), or not (Real)."
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_characteristic_axis_turns",  NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='characteristic_readout',
                              tooltip="Sets the relative phase of the sigma_x operation used to define the real/imag axis of the characteristic function.")

        # characteristic readout: configure readout grid
        self.setattr_argument("phases_pulse4_cat_turns",    PYONValue([0., 0.]), group='characteristic_readout', tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_pulse4_cat_phase",    EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'), group="characteristic_readout")
        self.setattr_argument("time_pulse4_cat_x_us_list",  Scannable(
                                                                default=[
                                                                    RangeScan(-50, 50, 10, randomize=True),
                                                                    ExplicitScan([50]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5), group="characteristic_readout")
        self.setattr_argument("time_pulse4_cat_y_us_list", Scannable(
                                                                default=[
                                                                    RangeScan(-50, 50, 10, randomize=True),
                                                                    ExplicitScan([50]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5), group="characteristic_readout")

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
        self._prepare_waveform()

        # # tmp remove
        # raise ValueError("Stop here")
        # # tmp remove

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
        self.att_displacement_mu = att_to_mu(self.att_displacement_db * dB)
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz
        self.freq_carrier_hz = self.freq_carrier_mhz * MHz
        self.freq_squeezing_base_hz = self.freq_squeezing_khz * kHz + self.freq_global_offset_hz

        # superresolution/characteristic reconstruction special
        self.time_gkp_stop_mu = self.core.seconds_to_mu(self.time_gjp_stop_us * us)
        # ensure pre-emptive stop time is multiple of phaser frame period
        self.time_gkp_stop_mu = ((self.time_gkp_stop_mu // self.phaser_eggs.t_frame_mu) *
                                             self.phaser_eggs.t_frame_mu)

    def _prepare_characteristic_reconstruction(self):
        """
        Prepare values for Characteristic Reconstruction.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS - DEFAULTS
        '''
        # defaults - singlepass AOM
        self.freq_singlepass_default_ftw_list = [self.singlepass0.frequency_to_ftw(freq_mhz * MHz)
                                                 for freq_mhz in self.freq_singlepass_default_mhz_list]
        self.ampl_singlepass_default_asf_list = [self.singlepass0.amplitude_to_asf(ampl_asf / 100.)
                                                 for ampl_asf in self.ampl_singlepass_default_pct_list]
        self.att_singlepass_default_mu_list =   [att_to_mu(att_db * dB)
                                                 for att_db in self.att_singlepass_default_db_list]

        # defaults - doublepass AOM
        self.ampl_doublepass_default_asf =  self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)
        self.att_doublepass_default_mu =    att_to_mu(self.att_doublepass_default_db * dB)

        # defaults - sigma_x pulses
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.att_sigmax_mu =    att_to_mu(self.att_sigmax_db * dB)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # defaults - cat
        self.freq_cat_center_ftw =  self.qubit.frequency_to_ftw(self.freq_cat_center_mhz * MHz)
        self.freq_cat_secular_ftw = self.qubit.frequency_to_ftw(self.freq_cat_secular_khz * kHz)
        self.ampls_cat_asf =    np.array([self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=np.int32)
        self.atts_cat_mu =      np.array([att_to_mu(att_db * dB)
                                          for att_db in self.atts_cat_db], dtype=np.int32)

        '''
        CONVERT VALUES TO MACHINE UNITS - PULSES
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

        # create sampling grid in radial coordinates
        vals_pulse4_mu_pow_list = np.array([
            [
                self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                self.singlepass0.turns_to_pow(np.arctan2(y_us, x_us) / (2. * np.pi))
            ]
            for x_us in self.time_pulse4_cat_x_us_list
            for y_us in self.time_pulse4_cat_y_us_list
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
            if np.sum(self.ampl_squeezing_frac) >= 100.:
                raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
            if np.abs(self.freq_squeezing_khz*kHz +  self.freq_global_offset_mhz * MHz) > 10 * MHz:
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_carrier_mhz * MHz)
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_carrier_mhz * MHz)
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

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
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.pulseshaper_vals = None        # store compiled waveforms from pulseshaper
        self.pulseshaper_id =   np.int32(0) # store waveform ID for pulseshaper

        # set up blocks for pulse sequence
        num_blocks = 1

        # set up the spin echo wizard generally
        # note: time_pulse_us is amount of time for each block
        self.spinecho_wizard.time_pulse_us =                self.squeezing_us / num_blocks
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           False
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        False
        self.spinecho_wizard.time_delay_spinecho_us =       250

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence & set amplitudes
        _sequence_blocks = np.zeros((num_blocks, 1, 2), dtype=float)
        _sequence_blocks[:, :, 0] = np.array(self.ampl_squeezing_frac, dtype=float)

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_ns = (self.core.mu_to_seconds(self.phaser_eggs.t_sample_mu) * ns) * np.arange(1)
        phase_osc_update_delay_turns_list = self.freq_squeezing_base_hz * t_update_delay_ns
        _sequence_blocks[:, :, 1] += np.array(self.phase_squeezing_turns) + phase_osc_update_delay_turns_list

        # create waveform
        self.spinecho_wizard.sequence_blocks = _sequence_blocks
        self.spinecho_wizard.calculate_pulseshape()
        self.spinecho_wizard.compile_waveform()

        # get waveform data and store in holding structure
        self.pulseshaper_vals = self.spinecho_wizard.get_waveform()

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
        # set up qubit beam for DMA sequences
        self.qubit.set_att_mu(self.att_doublepass_default_mu)
        delay_mu(10000)

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
            self.singlepass0.cpld.io_update.pulse_mu(8)
            delay_mu(8000)

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)

        '''INITIALIZATION'''
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
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
        delay_mu(1000000)
        # load waveform DMA handles
        self.pulse_shaper.waveform_load()
        delay_mu(250000)

        # used to check_termination more frequently
        _loop_iter = 0

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                time_char_read_mu =        config_vals[0]
                phase_char_read_pow =      np.int32(config_vals[1])
                characteristic_axis_bool =  bool(config_vals[2])
                self.core.break_realtime()

                # set phases for characteristic readout
                phase_char_read_pow_list = [
                    self.phases_pulse4_cat_pow[0] + self.phases_pulse4_cat_update_dir[0] * phase_char_read_pow,
                    self.phases_pulse4_cat_pow[1] + self.phases_pulse4_cat_update_dir[1] * phase_char_read_pow,
                ]

                # set phaser frequency
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    self.freq_carrier_hz - self.phaser_eggs.freq_center_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [self.freq_squeezing_hz, -self.freq_squeezing_hz, 0., 0., 0.],
                    self.phase_ch1_turns
                )
                self.core.break_realtime()

                '''MOTIONAL STATE PREPARATION'''
                # get current time
                t_phaser_start_mu = now_mu()

                # initialize ion in S-1/2 state & sideband cool to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # set target profile to ensure we run correctly (for characteristic reconstruction)
                self.qubit.set_profile(self.profile_729_target)
                self.qubit.cpld.io_update.pulse_mu(8)

                # apply squeezing
                if self.enable_phaser:
                    t_phaser_start_mu = self.phaser_run(self.pulseshaper_id)

                '''CHARACTERISTIC RECONSTRUCTION'''
                # prepare spin state for characteristic readout
                # note: need to set correct profile b/c might be stuck in SBC quench params
                self.pump.readout()
                self.repump_qubit.on()
                delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                self.repump_qubit.off()

                # sigma_x to select axis (does dummy if characteristic_axis_bool is False)
                self.pulse_sigmax(t_phaser_start_mu, self.phase_characteristic_axis_pow, characteristic_axis_bool)

                # bichromatic pulse + state detection
                self.pulse_bichromatic(t_phaser_start_mu, time_char_read_mu, phase_char_read_pow_list,
                                       self.freq_cat_center_ftw, self.freq_cat_secular_ftw)
                self.readout_subsequence.run_dma()

                '''LOOP CLEANUP'''
                # retrieve results & update dataset
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(characteristic_axis_bool, counts_res,
                                    time_char_read_mu, phase_char_read_pow)
                self.core.break_realtime()

                # resuscitate ion & detect deaths
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(counts_res)
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
            self.singlepass0.cpld.io_update.pulse_mu(8)
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

        ### CHARACTERISTIC SPECIAL ###
        if self.enable_cutoff:
            time_stop_mu = now_mu()
            # run a pre-emptive stop (atts + switches)
            at_mu(t_start_mu + self.time_superresolution_stop_mu)
            # self.phaser_eggs.channel[0].set_att_mu(0x00)
            # delay_mu(self.phaser_eggs.t_sample_mu)
            # self.phaser_eggs.channel[1].set_att_mu(0x00)
            self.phaser_eggs.ch0_amp_sw.off()
            self.phaser_eggs.ch1_amp_sw.off()

            at_mu(time_stop_mu)
        ### CHARACTERISTIC SPECIAL ###

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
        self.core.break_realtime()

        # record phaser pulse sequence and save returned waveform ID
        delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
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
        self.qubit.cpld.io_update.pulse_mu(8)

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_sigmax_mu << (0 * 8)) |
                (self.att_singlepass_default_mu_list[0] << (1 * 8)) |
                (self.att_singlepass_default_mu_list[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run sigmax pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        if is_real:
            self.qubit.on()
        else:
            self.qubit.off()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()

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
            freq_carrier_ftw, asf=self.ampl_doublepass_default_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
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

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_doublepass_default_mu << (0 * 8)) |
                (self.atts_cat_mu[0] << (1 * 8)) |
                (self.atts_cat_mu[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run bichromatic pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

