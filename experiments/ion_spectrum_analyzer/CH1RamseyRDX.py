from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS
from numpy import int32, int64, array, zeros, copy, arange, mean

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM,
    SidebandReadout, QubitRAP
)
from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper2 import PhaserPulseShaper2, PULSESHAPER_MAX_WAVEFORMS
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes


class CH1RamseyRDX(LAXExperiment, Experiment):
    """
    Experiment: CH1 Ramsey RDX

    Supports lots of easily configurable parameter scanning for phaser.
    """
    name = 'CH1 Ramsey RDX'
    kernel_invariants = {
        # hardware values - phaser
        'freq_sweep_hz_list', 'phase_sweep_turns_list', 'freq_osc_base_hz_list', 'freq_global_offset_hz',
        'waveform_index_to_phase_sweep_turns', 'phase_offsets',
        'waveform_index_to_pulseshaper_vals0', 'waveform_index_to_pulseshaper_vals1',

        # hardware values - readout
        'att_rap_mu', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu', 'freq_readout_ftw_list',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP',
        'spinecho_wizard', 'pulse_shaper',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',
        '_num_phaser_oscs',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=60, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config", BooleanValue(default=True))
        self.setattr_argument("sub_repetitions", NumberValue(default=1, precision=0, step=1, min=1, max=500),
                              tooltip="Loop over the same config for a given number of sub_repetitions."
                                      "Readout values are swept over for each sub_repetition.")
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="RAP"))

        self._num_phaser_oscs = 5   # number of phaser oscillators in use
        _argstr = "CH1Ram"  # string to use for argument grouping

        # allocate relevant beam profiles
        self.profile_729_sb_readout =   0
        self.profile_729_SBC =          1
        self.profile_729_RAP =          2

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_sb_readout)
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # build experiment arguments
        # extra argument for SBR
        self.setattr_argument("time_readout_us_list", Scannable(
                                                            default=[
                                                                ExplicitScan([26.11]),
                                                                RangeScan(0, 800, 150, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group='sideband_readout')
        self._build_arguments_RAP()
        self._build_arguments_freq_phase()
        self._build_arguments_waveform()

        # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def _build_arguments_RAP(self):
        """
        Set experiment arguments for RAP.
        """
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=100.7394, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=72., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=400., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

    def _build_arguments_freq_phase(self):
        """
        Set specific arguments for the phaser's frequency & phase.
        """
        _argstr = "CH1Ram"  # string to use for argument grouping

        # configurable freqs
        self.setattr_argument("freq_heating_carrier_mhz_list", Scannable(
                                                                default=[
                                                                    ExplicitScan([86.]),
                                                                    CenterScan(100, 0.002, 0.0005,
                                                                               randomize=True),
                                                                    # CenterScan(83.20175, 0.002, 0.0005, randomize=True),
                                                                ],
                                                                global_min=0.005, global_max=4800, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="{}.freq".format(_argstr))
        self.setattr_argument("freq_sweep_arr", PYONValue([1., 0., 0., 0., 0.]),
                              group="{}.freq".format(_argstr),
                              tooltip="Defines how oscillator freqs should be adjusted for each value in freq_superresolution_sweep_khz_list.\n"
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the freq value, and osc_1 by -1x the freq value, with the rest untouched.\n"
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("freq_sweep_khz_list",    Scannable(
                                                            default=[
                                                                ExplicitScan([0.]),
                                                                CenterScan(0., 4., 0.1, randomize=True),
                                                            ],
                                                            global_min=-10000, global_max=10000, global_step=10,
                                                            unit="kHz", scale=1, precision=6
                                                        ), group = "{}.freq".format(_argstr))

        # configurable phases
        self.setattr_argument("target_phase_sweep", EnumerationValue(['ch0+ch1', 'ch1'], default='ch0+ch1'),
                              group="{}.phase".format(_argstr),
                              tooltip="Choose whether phase sweep is applied uniformly to oscillators on both channels, "
                                      "or only oscillators on CH1.\n"
                                      "Phase sweep is applied to BOTH Ramsey stages.")
        self.setattr_argument("phase_sweep_arr", PYONValue([0., 0., 0., 0., 0.]), group="{}.phase".format(_argstr),
                              tooltip="Defines how oscillator phases should be adjusted for each value in phase_sweep_turns_list.\n"
                                      "Depending on value of target_phase_sweep, this will be applied to both CH0 and CH1, or only CH1.\n"
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("phase_sweep_turns_list", Scannable(
                                                            default=[
                                                                ExplicitScan([0.]),
                                                                RangeScan(0, 1.0, 26, randomize=True),
                                                                ExplicitScan([0.]),
                                                            ],
                                                            global_min=0.0, global_max=1.0, global_step=1,
                                                            unit="turns", scale=1, precision=5
                                                        ), group = "{}.phase".format(_argstr))
        self.setattr_argument("phase_global_ch1_duc_turns", NumberValue(0., min=0.0, max=1.0, step=1,
                                                                        unit="turns", scale=1, precision=5),
                                                                        group="{}.phase".format(_argstr))

    def _build_arguments_waveform(self):
        """
        Set specific arguments for the phaser waveform.
        """
        _argstr = "CH1Ram"  # string to use for argument grouping

        # RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.shape'.format(_argstr))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.shape'.format(_argstr))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.shape'.format(_argstr))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.shape'.format(_argstr))

        # custom waveform specification - general
        self.setattr_argument("time_heating_us",   NumberValue(default=1e3, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.wav".format(_argstr))
        self.setattr_argument("att_phaser_db_list", Scannable(
                                                        default=[
                                                            ExplicitScan([31.5]),
                                                            RangeScan(15., 31.5, 33, randomize=True),
                                                        ],
                                                        global_min=0.0, global_max=31.5, global_step=0.5,
                                                        unit="dB", scale=1, precision=1
                                                    ), group = "{}.wav".format(_argstr))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.wav".format(_argstr),
                              tooltip="Attempt to move NCO/TRF leakage outside of target band by shifting DUC and oscillators to compensate.")
        self.setattr_argument("freq_osc_khz_list",      PYONValue([702.581, 0., 0., 0., 0.]), group="{}.wav".format(_argstr))
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0., 0.]), group="{}.wav".format(_argstr),
                              tooltip="Relative phases between each oscillator. Applies equally to CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0., 0., 0.]), group="{}.wav".format(_argstr),
                              tooltip="Individual CH1 offsets for each oscillator. Obviously, applies only to CH1.\n"
                                      "This is in addition to the CH1 global offset, as well as any CH1 sweeps.")

        # custom waveform specification - CH1 Ramsey-specific
        self.setattr_argument('ampl_ch0_stage_0',   PYONValue([40., 40., 0., 0., 0.]), group= "{}.wav".format(_argstr),
                              tooltip="Amplitudes for CH0 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch0_stage_1',   PYONValue([40., 40., 0., 0., 0.]), group="{}.wav".format(_argstr),
                              tooltip="Amplitudes for CH0 oscillators during second stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_0',   PYONValue([40., 40., 0., 0., 0.]), group= "{}.wav".format(_argstr),
                              tooltip="Amplitudes for CH1 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_1',   PYONValue([40., 40., 5, 0., 0.]), group="{}.wav".format(_argstr),
                              tooltip="Amplitudes for CH1 oscillators during second stage of Ramsey.")


    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''GENERAL SETUP'''
        self._prepare_argument_checks()

        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_oscillators_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper2(self, array(self.phase_osc_ch1_offset_turns))


        '''SUBSEQUENCE PARAMETERS'''
        # prepare RAP arguments
        self.att_rap_mu = att_to_mu(self.att_rap_db * dB)
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        # configure readout method
        if self.readout_type == 'RAP':
            self.enable_RAP = True
            self.freq_readout_ftw_list = array([self.freq_rap_center_ftw], dtype=int32)
            time_readout_mu_list = [self.time_rap_mu]
        elif self.readout_type == 'Sideband Ratio':
            self.enable_RAP = False
            self.freq_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
            time_readout_mu_list = [self.core.seconds_to_mu(time_us * us) for time_us in self.time_readout_us_list]
        else:
            raise ValueError("Invalid readout type. Must be one of (Sideband Ratio, RAP).")


        '''HARDWARE VALUES - CONFIG'''
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz
        # ensure att_phaser_mu_list rounded to 0.5dB (min att step size) to prevent spamming
        att_phaser_mu_list = set(att_to_mu(round(att_db * 2) / 2 * dB) for att_db in self.att_phaser_db_list)

        # convert build arguments to appropriate values and format as numpy arrays
        freq_carrier_hz_list = array(list(self.freq_heating_carrier_mhz_list)) * MHz
        self.freq_sweep_hz_list = array(list(self.freq_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_hz
        self.phase_sweep_turns_list = list(self.phase_sweep_turns_list)
        self.phase_osc_ch1_offset_turns = array(self.phase_osc_ch1_offset_turns)
        self.freq_sweep_arr = array(self.freq_sweep_arr, dtype=float)
        self.phase_sweep_arr = array(self.phase_sweep_arr, dtype=float)


        '''CONFIGURE SWEEP BEHAVIOR'''
        # workaround: pre-declare phase_offsets b/c stack memory issue (https://github.com/m-labs/artiq/issues/1520)
        if self.target_phase_sweep == "ch1":
            self.phase_offsets = array([self.phase_osc_ch1_offset_turns + self.phase_sweep_arr * phase_val_turns
                                           for phase_val_turns in self.phase_sweep_turns_list])
        elif self.target_phase_sweep == "ch0+ch1":
            self.phase_offsets = array([self.phase_osc_ch1_offset_turns])


        '''CREATE EXPERIMENT CONFIG'''
        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_sweep_turns = arange(len(self.phase_sweep_turns_list))

        # create config data structure
        self.config_experiment_list = create_experiment_config(
            freq_carrier_hz_list,
            self.freq_sweep_hz_list, self.waveform_index_to_phase_sweep_turns,
            time_readout_mu_list, att_phaser_mu_list,
            shuffle_config=self.randomize_config, config_type=int32
        )

        self._prepare_waveform() # configure waveform via pulse shaper & spin echo wizard

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
            elif any((sum(ampl_arg) >= 100. for ampl_arg in ampl_arg_tuple)):
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

        max_osc_freq_hz = (max(list(self.freq_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        min_osc_freq_hz = (min(list(self.freq_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = array(list(self.freq_heating_carrier_mhz_list)) * MHz
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
        if not isinstance(self.phase_sweep_arr, list):
            raise ValueError("Error: phaser oscillator on/off must be a list.")
        elif len(self.phase_sweep_arr) != self._num_phaser_oscs:
            raise ValueError("Error: only {:d} oscillators to change phase.".format(self._num_phaser_oscs))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for pulse waveforms
        self.waveform_index_to_pulseshaper_vals0 =  list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_vals1 =  list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_id =     zeros(len(self.phase_sweep_turns_list),
                                                             dtype=int32) # store pulseshaper waveform ID
        num_blocks = 2  # set up blocks for pulse sequence

        '''DESIGN WAVEFORM SEQUENCE'''
        # create separate bare waveform block sequences for CH0 and CH1
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        _osc_vals_ch0 = zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_ch1 = zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)

        # update CH0 and CH1 sequences w/relevant amplitudes
        _osc_vals_ch0[:, :, 0] = array([self.ampl_ch0_stage_0, self.ampl_ch0_stage_1])
        _osc_vals_ch1[:, :, 0] = array([self.ampl_ch1_stage_0, self.ampl_ch1_stage_1])

        # phase track oscillator updates to account for 40ns sample period
        t_update_delay_s_list = array([0, 40e-9, 80e-9, 80e-9, 120e-9])
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_sweep_arr * mean(self.freq_sweep_hz_list)) *
                t_update_delay_s_list
        )
        _osc_vals_ch0[:, :, 1] += array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list
        _osc_vals_ch1[:, :, 1] += array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list


        '''COMPILE WAVEFORM SEQUENCE'''
        for phase in self.phase_sweep_turns_list:
            # create local copy of _sequence_blocks
            # note: no need to deep copy b/c it's filled w/immutables
            # note: have to obtain different copies so they don't point to same object and overwrite it
            _osc_vals_local_ch0 = copy(_osc_vals_ch0)
            _osc_vals_local_ch1 = copy(_osc_vals_ch1)

            # apply oscillator phase sweep
            if self.target_phase_sweep == "ch0+ch1":
                _osc_vals_local_ch0[:, :, 1] += self.phase_sweep_arr * phase
                _osc_vals_local_ch1[:, :, 1] += self.phase_sweep_arr * phase

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
        return (self.repetitions * self.sub_repetitions *
                len(self.config_experiment_list) * len(self.freq_readout_ftw_list),
                7)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # set beams to rescue while we wait (long initialize for phaser-type exps)
        self.initialize_subsequence.slack_rescue()

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
        self.pulse_shaper.waveform_load() # load waveform DMA handles
        _loop_iter = 0 # used to check_termination more frequently

        '''MAIN LOOP'''
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                carrier_freq_hz =   config_vals[0]
                freq_sweep_hz =     config_vals[1]
                phase_sweep_idx =   int32(config_vals[2])
                time_readout_mu =   int64(config_vals[3])
                att_phaser_mu =     config_vals[4]

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


                '''SUB-REPETITION IMPLEMENTATION'''
                for _subrep_num in range(self.sub_repetitions):
                    # sub-repetitions requires we run a bunch of sub_reps at the same config
                    # however, we still want to sweep the readout parameters
                    for freq_readout_ftw in self.freq_readout_ftw_list:

                        # configure readout
                        self.core.break_realtime()
                        if not self.enable_RAP:
                            self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                              profile=self.profile_729_sb_readout,
                                              phase_mode=PHASE_MODE_CONTINUOUS)
                            delay_mu(25000)


                        '''STATE PREPARATION'''
                        self.initialize_subsequence.run_dma()
                        self.sidebandcool_subsequence.run_dma()

                        '''HEATING'''
                        self.phaser_run(waveform_id, att_phaser_mu)

                        '''READOUT'''
                        if self.enable_RAP:
                            self.qubit.set_att_mu(self.att_rap_mu)
                            self.rap_subsequence.run_rap(time_readout_mu)
                        else:
                            self.sidebandreadout_subsequence.run_time(time_readout_mu)
                        self.readout_subsequence.run_dma()

                        '''CLEANUP'''
                        # clean up
                        self.rescue_subsequence.resuscitate()
                        self.initialize_subsequence.slack_rescue()

                        # retrieve results & update dataset
                        counts = self.readout_subsequence.fetch_count()
                        self.rescue_subsequence.detect_death(counts)
                        self.update_results(freq_readout_ftw,
                                            counts,
                                            carrier_freq_hz,
                                            freq_sweep_hz,
                                            phase_sweep_turns,
                                            time_readout_mu,
                                            att_phaser_mu)

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
    def phaser_run(self, waveform_id: TInt32, att_mu: TInt32) -> TNone:
        """
        Run the main  pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        """
        # START/SETUP
        self.phaser_eggs.phaser_setup(att_mu, att_mu)

        # RUN
        self.phaser_eggs.reset_duc_phase() # reset DUC for deterministic start phase
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
            # get waveforms for given parameters
            # note: use sync RPC to reduce significant overhead of direct data transfer
            _wav_data_ampl0, _wav_data_phas0, _wav_data_time0 = self._get_compiled_waveform_ch0(i)
            _wav_data_ampl1, _wav_data_phas1, _wav_data_time1 = self._get_compiled_waveform_ch1(i)
            if self.target_phase_sweep == 'ch1':
                self.pulse_shaper.phase_offsets_turns = self.phase_offsets[i]

            # record phaser pulse sequence and save returned waveform ID
            # note: no need to add slack b/c waveform_record does it for us
            self.waveform_index_to_pulseshaper_id[i] = self.pulse_shaper.waveform_record(
                _wav_data_ampl0, _wav_data_ampl1,
                _wav_data_phas0, _wav_data_phas1,
                _wav_data_time0
            )

    @rpc
    def _get_compiled_waveform_ch0(self, wav_idx: TInt32) -> TTuple([TArray(TFloat, 2),
                                                                     TArray(TFloat, 2),
                                                                     TArray(TInt64, 1)]):
        """
        Return compiled waveform values - CH0.
        By returning the large waveform arrays via RPC, we avoid all-at-once large data transfers,
            speeding up experiment compilation and transfer to Kasli.
        :param wav_idx: the index of the compiled waveform to retrieve.
        """
        return self.waveform_index_to_pulseshaper_vals0[wav_idx]

    @rpc
    def _get_compiled_waveform_ch1(self, wav_idx: TInt32) -> TTuple([TArray(TFloat, 2),
                                                                     TArray(TFloat, 2),
                                                                     TArray(TInt64, 1)]):
        """
        Return compiled waveform values - CH1.
        By returning the large waveform arrays via RPC, we avoid all-at-once large data transfers,
            speeding up experiment compilation and transfer to kasli.
        :param wav_idx: the index of the compiled waveform to retrieve.
        """
        return self.waveform_index_to_pulseshaper_vals1[wav_idx]


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

