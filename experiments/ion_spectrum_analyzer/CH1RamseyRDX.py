import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM,
    SidebandReadout, QubitRAP
)
from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper2 import PhaserPulseShaper2, PULSESHAPER_MAX_WAVEFORMS


class CH1RamseyRDX(LAXExperiment, Experiment):
    """
    Experiment: CH1 Ramsey RDX

    Supports lots of easily configurable parameter scanning for phaser.
    """
    name = 'CH1 Ramsey RDX'
    kernel_invariants = {
        # hardware values - phaser
        'att_heating_mu', 'freq_sweep_hz_list', 'phase_sweep_turns_list', 'freq_osc_base_hz_list',
        'waveform_index_to_pulseshaper_vals0', 'waveform_index_to_pulseshaper_vals1', 'freq_global_offset_hz',
        'waveform_index_to_phase_sweep_turns', 'phase_offsets',

        # hardware values - readout
        'att_rap_mu', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP',
        'spinecho_wizard', 'pulse_shaper',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',
        '_num_phaser_oscs'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=1000, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="Sideband Ratio"))

        self._num_phaser_oscs = 5   # number of phaser oscillators in use
        _argstr = "CH1Ram"  # string to use for argument grouping

        # allocate relevant beam profiles
        self.profile_729_sb_readout =   0
        self.profile_729_SBC =          1
        self.profile_729_RAP =          2

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )

        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_sb_readout)
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # Sideband Readout - extra argument
        self.setattr_argument("time_readout_us_list", Scannable(
                                                            default=[
                                                                ExplicitScan([26.11]),
                                                                RangeScan(0, 800, 150, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group='sideband_readout')

        # RAP-based readout
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=100.7367, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

        # configurable freq & sweeps
        self.setattr_argument("freq_heating_carrier_mhz_list", Scannable(
                                                                default=[
                                                                    ExplicitScan([70.]),
                                                                    CenterScan(100, 0.002, 0.0005,
                                                                               randomize=True),
                                                                    # CenterScan(83.20175, 0.002, 0.0005, randomize=True),
                                                                ],
                                                                global_min=0.005, global_max=4800, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("freq_sweep_arr", PYONValue([-1., 1., 0., 0., 0.]), group="{}.freq_phase_sweep".format(_argstr),
                              tooltip="Defines how oscillator freqs should be adjusted for each value in freq_superresolution_sweep_khz_list."
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the freq value, and osc_1 by -1x the freq value, with the rest untouched."
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("freq_sweep_khz_list",    Scannable(
                                                            default=[
                                                                ExplicitScan([0.]),
                                                                CenterScan(0., 4., 0.1, randomize=True),
                                                            ],
                                                            global_min=-10000, global_max=10000, global_step=10,
                                                            unit="kHz", scale=1, precision=6
                                                        ), group = "{}.freq_phase_sweep".format(_argstr))

        self.setattr_argument("target_phase_sweep", EnumerationValue(['osc0', 'osc1', 'osc2', 'osc3', 'osc4', 'ch1'], default='osc0'),
                              group = "{}.freq_phase_sweep".format(_argstr),
                              tooltip="Phase sweep is applied to BOTH Ramsey stages.")
        self.setattr_argument("phase_sweep_turns_list", Scannable(
                                                            default=[
                                                                ExplicitScan([0.]),
                                                                RangeScan(0, 1.0, 26, randomize=True),
                                                                ExplicitScan([0.]),
                                                            ],
                                                            global_min=0.0, global_max=1.0, global_step=1,
                                                            unit="turns", scale=1, precision=3
                                                        ), group = "{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("phase_global_ch1_duc_turns", NumberValue(0., min=0.0, max=1.0, step=1,
                                                                        unit="turns", scale=1, precision=3),
                                                                        group="{}.freq_phase_sweep".format(_argstr))

        # RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        # custom waveform specification - general
        self.setattr_argument("time_heating_us",   NumberValue(default=1e3, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.waveform".format(_argstr))
        self.setattr_argument("att_heating_db",    NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format(_argstr))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Attempt to move NCO/TRF leakage outside of target band by shifting DUC and oscillators to compensate.")
        self.setattr_argument("freq_osc_khz_list",      PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format(_argstr),
                              tooltip="Relative phases between each oscillator. Applies equally to CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format(_argstr),
                              tooltip="Individual CH1 offsets for each oscillator. Obviously, applies only to CH1."
                                      "This is in addition to the CH1 global offset, as well as any CH1 sweeps.")

        # custom waveform specification - CH1 Ramsey-specific
        self.setattr_argument('ampl_ch0_stage_0',   PYONValue([10., 10., 5., 0., 0.]), group= "{}.waveform".format(_argstr),
                              tooltip="Amplitudes for CH0 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch0_stage_1',   PYONValue([10., 10., 5., 0., 0.]), group="{}.waveform".format(_argstr),
                              tooltip="Amplitudes for CH0 oscillators during second stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_0',   PYONValue([10., 10., 5., 0., 0.]), group= "{}.waveform".format(_argstr),
                              tooltip="Amplitudes for CH1 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_1',   PYONValue([10., 10., 5, 0., 0.]), group="{}.waveform".format(_argstr),
                              tooltip="Amplitudes for CH1 oscillators during second stage of Ramsey.")
        self.setattr_argument('osc_ch1_sweep', PYONValue([0., 0., 0., 0., 0.]), group= "{}.waveform".format(_argstr),
                              tooltip="Sweep array when target_phase_sweep is set to 'ch1'. Each oscillator ")

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE & VALIDATE INPUTS'''
        self._prepare_argument_checks()

        '''SUBSEQUENCE PARAMETERS'''
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_oscillators_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper2(self, np.array(self.phase_osc_ch1_offset_turns))

        # prepare RAP arguments
        self.att_rap_mu = att_to_mu(self.att_rap_db * dB)
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        # configure readout method
        if self.readout_type == 'RAP':
            self.enable_RAP = True
            freq_sideband_readout_ftw_list =    np.array([self.freq_rap_center_ftw], dtype=np.int32)
            time_readout_mu_list =              np.array([self.time_rap_mu], dtype=np.int64)
        elif self.readout_type == 'Sideband Ratio':
            self.enable_RAP = False
            freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
            time_readout_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                             for time_us in self.time_readout_us_list])
        else:
            raise ValueError("Invalid readout type. Must be one of (Sideband Ratio, RAP).")

        '''HARDWARE VALUES - CONFIG'''
        # convert hardware values to convenient units
        self.att_heating_mu = att_to_mu(self.att_heating_db * dB)
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # convert build arguments to appropriate values and format as numpy arrays
        freq_carrier_hz_list = np.array(list(self.freq_heating_carrier_mhz_list)) * MHz
        self.freq_sweep_hz_list = np.array(list(self.freq_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = np.array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_hz
        self.phase_sweep_turns_list = np.array(list(self.phase_sweep_turns_list))
        self.phase_osc_ch1_offset_turns = np.array(self.phase_osc_ch1_offset_turns)
        self.osc_ch1_sweep = np.array(self.osc_ch1_sweep)
        self.freq_sweep_arr = np.array(self.freq_sweep_arr, dtype=float)


        '''CONFIGURE SWEEP BEHAVIOR'''
        # create pre-declared phase_offsets list as workaround for artiq stack memory issue
        # see: https://github.com/m-labs/artiq/issues/1520
        self.phase_offsets = np.zeros((len(self.phase_sweep_turns_list), self._num_phaser_oscs))
        for i, phase_val_turns in enumerate(self.phase_sweep_turns_list):
            self.phase_offsets[i, :] = self.phase_osc_ch1_offset_turns + phase_val_turns * self.osc_ch1_sweep

        '''CREATE EXPERIMENT CONFIG'''
        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_sweep_turns = np.arange(len(self.phase_sweep_turns_list))

        # create config data structure
        self.config_experiment_list = np.zeros((
            len(freq_sideband_readout_ftw_list) *
            len(freq_carrier_hz_list) *
            len(self.freq_sweep_hz_list) *
            len(self.phase_sweep_turns_list) *
            len(time_readout_mu_list),
        5), dtype=float)
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_experiment_list[:, [1, 2, -3, -1, 0]] = np.stack(np.meshgrid(
            freq_carrier_hz_list,
            self.freq_sweep_hz_list,
            self.waveform_index_to_phase_sweep_turns,
            time_readout_mu_list,
            freq_sideband_readout_ftw_list),
        -1).reshape(-1, 5)
        # randomize_config always enabled lol
        np.random.shuffle(self.config_experiment_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # check that input amplitude/phase arrays are valid
        ampl_arg_tuple = (self.ampl_ch0_stage_0, self.ampl_ch0_stage_1, self.ampl_ch1_stage_0, self.ampl_ch1_stage_0)
        ampls_are_list = all((isinstance(ampl_arg, list) for ampl_arg in ampl_arg_tuple))
        if ampls_are_list:
            if any((len(ampl_arg) != self._num_phaser_oscs for ampl_arg in ampl_arg_tuple)):
                raise ValueError("Error: phaser oscillator amplitude arrays must have length {:d}.".format(self._num_phaser_oscs))
            elif any((np.sum(ampl_arg) >= 100. for ampl_arg in ampl_arg_tuple)):
                raise ValueError("Error: phaser oscillator amplitudes must sum < 100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude arrays must be lists.")

        if isinstance(self.phase_osc_turns_list, list):
            if len(self.phase_osc_turns_list) != self._num_phaser_oscs:
                raise ValueError("Error: phaser oscillator phase array must have length {:d}.".format(self._num_phaser_oscs))
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        if ((not isinstance(self.freq_osc_khz_list, list)) or
                (len(self.freq_osc_khz_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator frequency array must be list of length {:d}.".format(self._num_phaser_oscs))

        max_osc_freq_hz = (np.max(list(self.freq_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        min_osc_freq_hz = (np.min(list(self.freq_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_heating_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        # ensure that sweep targets are valid
        if not (isinstance(self.freq_sweep_arr, list) and (len(self.freq_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid freq_sweep_arr: {:}."
                             "freq_sweep_arr must be list of length {:d}.".format(self.freq_sweep_arr,
                                                                                      self._num_phaser_oscs))

        # ensure special CH1 ramsey stuff is valid
        if not isinstance(self.osc_ch1_sweep, list):
            raise ValueError("Error: phaser oscillator on/off must be a list.")
        elif len(self.osc_ch1_sweep) != self.self._num_phaser_oscs:
            raise ValueError("Error: only {:d} oscillators to change phase.".format(self._num_phaser_oscs))
        elif np.max(self.osc_ch1_sweep) > 1:
            raise ValueError("Error: can only turn oscillator on (1)/off (0)")
        elif np.max(self.osc_ch1_sweep) < 0:
            raise ValueError("Error: can only turn oscillator on (1)/off (0)")

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for pulse waveforms
        self.waveform_index_to_pulseshaper_vals0 =  list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_vals1 =  list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_id =     np.zeros(len(self.phase_sweep_turns_list),
                                                             dtype=np.int32) # store pulseshaper waveform ID

        num_blocks = 2  # set up blocks for pulse sequence

        # prepare phase sweep
        if self.target_phase_sweep == "osc0":
            phas_update_arr = np.array([1., 0., 0., 0., 0.])
        elif self.target_phase_sweep == "osc1":
            phas_update_arr = np.array([0., 1., 0., 0., 0.])
        elif self.target_phase_sweep == "osc2":
            phas_update_arr = np.array([0., 0., 1., 0., 0.])
        elif self.target_phase_sweep == "osc3":
            phas_update_arr = np.array([0., 0., 0., 1., 0.])
        elif self.target_phase_sweep == "osc4":
            phas_update_arr = np.array([0., 0., 0., 0., 1.])
        elif self.target_phase_sweep == "ch1":
            phas_update_arr = np.array([0., 0., 0., 0., 0.])
        else:
            raise ValueError("Invalid phase sweep type.")


        '''DESIGN WAVEFORM SEQUENCE'''
        # create separate bare waveform block sequences for CH0 and CH1
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        _osc_vals_ch0 = np.zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_ch1 = np.zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)

        # update CH0 and CH1 sequences w/relevant amplitudes
        _osc_vals_ch0[:, :, 0] = np.array([self.ampl_ch0_stage_0, self.ampl_ch0_stage_1])
        _osc_vals_ch1[:, :, 0] = np.array([self.ampl_ch1_stage_0, self.ampl_ch1_stage_1])

        # phase track oscillator updates to account for 40ns sample period
        t_update_delay_s_list = np.array([0, 40e-9, 80e-9, 80e-9, 120e-9])
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_sweep_arr * np.mean(self.freq_sweep_hz_list)) *
                t_update_delay_s_list
        )
        _osc_vals_ch0[:, :, 1] += np.array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list
        _osc_vals_ch1[:, :, 1] += np.array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list


        '''COMPILE WAVEFORM SEQUENCE'''
        for i, phase in enumerate(self.phase_sweep_turns_list):
            # create local copy of _sequence_blocks
            # note: no need to deep copy b/c it's filled w/immutables
            # note: have to obtain different copies so they don't point to same object and overwrite it
            _osc_vals_local_ch0 = np.copy(_osc_vals_ch0)
            _osc_vals_local_ch1 = np.copy(_osc_vals_ch1)
            _osc_vals_local_ch0[:, :, 1] += phas_update_arr * phase
            _osc_vals_local_ch1[:, :, 1] += phas_update_arr * phase

            _sequence_blocks_local_ch0 = [
                {
                    "oscillator_parameters": _osc_vals_local_ch0[i],
                    "config": {
                        "time_us": self.time_heating_us / num_blocks,
                        "pulse_shaping": self.enable_pulse_shaping,
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping,
                            "pulse_shape_falling": self.enable_pulse_shaping,
                            "sample_rate_khz": self.freq_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_pulse_shape_rolloff_us
                        }
                    }
                } for i in range(num_blocks)
            ]
            _sequence_blocks_local_ch1 = [
                {
                    "oscillator_parameters": _osc_vals_local_ch1[i],
                    "config": {
                        "time_us": self.time_heating_us / num_blocks,
                        "pulse_shaping": self.enable_pulse_shaping,
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping,
                            "pulse_shape_falling": self.enable_pulse_shaping,
                            "sample_rate_khz": self.freq_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_pulse_shape_rolloff_us
                        }
                    }
                } for i in range(num_blocks)
            ]

            # compile waveform and store in holding structure
            self.waveform_index_to_pulseshaper_vals0.append(self.spinecho_wizard.compile_waveform(_sequence_blocks_local_ch0))
            self.waveform_index_to_pulseshaper_vals1.append(self.spinecho_wizard.compile_waveform(_sequence_blocks_local_ch1))

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                6)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # configure RAP pulse
        if self.enable_RAP:
            self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
            delay_mu(50000)

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
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                carrier_freq_hz =   config_vals[1]
                freq_sweep_hz =     config_vals[2]
                phase_sweep_idx =   np.int32(config_vals[3])
                time_readout_mu =   np.int64(config_vals[4])

                # get corresponding phase and waveform ID from the index
                phase_sweep_turns = self.phase_sweep_turns_list[phase_sweep_idx]
                waveform_id = self.waveform_index_to_pulseshaper_id[phase_sweep_idx]

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_osc_base_hz_list + (freq_sweep_hz * self.freq_sweep_arr)
                self.core.break_realtime()
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    carrier_freq_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], freq_update_list[4]],
                    self.phase_global_ch1_duc_turns
                )

                # set qubit readout frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_sb_readout, phase_mode=PHASE_MODE_CONTINUOUS)
                self.core.break_realtime()


                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                '''HEATING'''
                self.phaser_run(waveform_id)

                '''READOUT'''
                if self.enable_RAP:
                    self.qubit.set_att_mu(self.att_rap_mu)
                    self.rap_subsequence.run_rap(time_readout_mu)
                else:
                    self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
                self.update_results(
                    freq_readout_ftw,
                    counts,
                    carrier_freq_hz,
                    freq_sweep_hz,
                    phase_sweep_turns,
                    time_readout_mu
                )

                '''LOOP CLEANUP'''
                # resuscitate ion & detect deaths
                self.core.break_realtime()
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(counts)

                # check termination more frequently in case reps are low
                if _loop_iter % 50 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed and support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    HELPER FUNCTIONS - PHASER
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main  pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        """
        # START/SETUP
        self.phaser_eggs.phaser_setup(self.att_heating_mu, self.att_heating_mu)

        # RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # STOP
        # stop all output & clean up hardware (e.g. amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each RSB phase
        for i in range(len(self.phase_sweep_turns_list)):
            # get waveform
            _wav_data_ampl0, _wav_data_phas0, _wav_data_time0 = self.waveform_index_to_pulseshaper_vals0[i]
            _wav_data_ampl1, _wav_data_phas1, _wav_data_time1 = self.waveform_index_to_pulseshaper_vals1[i]
            if self.target_phase_sweep == 'ch1':
                self.pulse_shaper.phase_offsets_turns = self.phase_offsets[i]
            self.core.break_realtime()

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            _wav_idx = self.pulse_shaper.waveform_record(_wav_data_ampl0, _wav_data_ampl1,
                                                         _wav_data_phas0, _wav_data_phas1,
                                                         _wav_data_time0)
            self.waveform_index_to_pulseshaper_id[i] = _wav_idx
            self.core.break_realtime()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

