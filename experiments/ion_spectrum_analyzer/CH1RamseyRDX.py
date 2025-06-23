import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM,
    SidebandReadout, QubitRAP
)
from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper2 import PhaserPulseShaper2

from sipyco import pyon


class CH1RamseyRDX(LAXExperiment, Experiment):
    """
    Experiment: CH1 Ramsey RDX

    Supports lots of easily configurable parameter scanning for phaser.
    """
    name = 'CH1 Ramsey RDX'
    kernel_invariants = {
        # hardware values - phaser
        'att_heating_mu', 'freq_sweep_hz_list', 'freq_update_arr', 'phase_sweep_turns_list', 'freq_osc_base_hz_list',
        'waveform_index_to_pulseshaper_vals0', 'waveform_index_to_pulseshaper_vals1',
        'waveform_index_to_phase_sweep_turns', 'phase_offsets',

        # hardware values - readout
        'att_rap_mu', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',

        # subharmonic specials
        'freq_global_offset_hz'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=1000, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="Sideband Ratio"))
        self.setattr_argument("randomize_config", BooleanValue(default=True))

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
                                                                ExplicitScan([50.89]),
                                                                RangeScan(0, 1500, 100, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group='sideband_readout')

        # RAP-based readout
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=101.0991, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=150., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=1000., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

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
                                                                    ), group="{}.freq_phase_sweep".format("SDR"))

        self.setattr_argument("target_freq_sweep",  EnumerationValue(['Secular', 'Probe', 'osc0', 'osc1', 'osc2', 'osc3', 'osc4'],
                                                                     default='Secular'), group = "{}.freq_phase_sweep".format("SDR"))
        self.setattr_argument("freq_sweep_khz_list",    Scannable(
                                                                            default=[
                                                                                # ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                                ExplicitScan([0.]),
                                                                                CenterScan(0., 4., 0.1, randomize=True),
                                                                            ],
                                                                            global_min=-10000, global_max=10000, global_step=10,
                                                                            unit="kHz", scale=1, precision=6
                                                                        ), group = "{}.freq_phase_sweep".format("SDR"))

        self.setattr_argument("target_phase_sweep", EnumerationValue(['osc0', 'osc1', 'osc2', 'osc3', 'osc4', 'ch1'], default='osc0'),
                              group = "{}.freq_phase_sweep".format("SDR"),
                              tooltip="Phase sweep is applied to BOTH Ramsey stages.")
        self.setattr_argument("phase_sweep_turns_list", Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        RangeScan(0, 1.0, 3, randomize=True),
                                                                        ExplicitScan([0.]),
                                                                    ],
                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                    unit="turns", scale=1, precision=3
                                                                ), group = "{}.freq_phase_sweep".format("SDR"))
        self.setattr_argument("phase_heating_global_ch1_turns_list", NumberValue(0.367, min=0.0, max=1.0, step=1,
                                                                        unit="turns", scale=1, precision=3),
                                                                        group="{}.freq_phase_sweep".format("SDR"))

        # RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format("SDR"))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format("SDR"))

        # custom waveform specification - general
        self.setattr_argument("time_heating_us",   NumberValue(default=1e3, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("att_heating_db",    NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format("SDR"))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.waveform".format("SDR"),
                              tooltip="Attempt to move NCO/TRF leakage outside of target band by shifting DUC and oscillators to compensate.")
        self.setattr_argument("freq_osc_khz_list",      PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format("SDR"))
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format("SDR"),
                              tooltip="Relative phases between each oscillator. Applies equally to CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0., 0., 0.]), group="{}.waveform".format("SDR"),
                              tooltip="Individual CH1 offsets for each oscillator. Obviously, applies only to CH1."
                                      "This is in addition to the CH1 global offset, as well as any CH1 sweeps.")

        # custom waveform specification - CH1 Ramsey-specific
        self.setattr_argument('ampl_ch0_stage_0',   PYONValue([10., 10., 5., 0., 0.]), group= "{}.waveform".format("SDR"),
                              tooltip="Amplitudes for CH0 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch0_stage_1',   PYONValue([10., 10., 5., 0., 0.]), group="{}.waveform".format("SDR"),
                              tooltip="Amplitudes for CH0 oscillators during second stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_0',   PYONValue([10., 10., 5., 0., 0.]), group= "{}.waveform".format("SDR"),
                              tooltip="Amplitudes for CH1 oscillators during first stage of Ramsey.")
        self.setattr_argument('ampl_ch1_stage_1',   PYONValue([10., 10., 5, 0., 0.]), group="{}.waveform".format("SDR"),
                              tooltip="Amplitudes for CH1 oscillators during second stage of Ramsey.")
        self.setattr_argument('osc_ch1_sweep', PYONValue([0., 0., 0., 0., 0.]), group= "{}.waveform".format("SDR"),
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
        self.freq_rap_center_ftw = hz_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = hz_to_ftw(self.freq_rap_dev_khz * kHz)
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
        # convert attenuation from dB to machine units
        self.att_heating_mu = att_to_mu(self.att_heating_db * dB)

        # convert freqs to Hz
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # convert build arguments to appropriate values and format as numpy arrays
        freq_carrier_hz_list = np.array(list(self.freq_heating_carrier_mhz_list)) * MHz
        self.freq_sweep_hz_list = np.array(list(self.freq_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = np.array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_hz
        self.phase_sweep_turns_list = np.array(list(self.phase_sweep_turns_list))
        self.phase_heating_global_ch1_turns_list = np.array(list([self.phase_heating_global_ch1_turns_list]))
        self.phase_osc_ch1_offset_turns = np.array(self.phase_osc_ch1_offset_turns)
        self.osc_ch1_sweep = np.array(self.osc_ch1_sweep)

        '''CONFIGURE SWEEP BEHAVIOR'''
        # implement variable freq sweep target
        # adjust oscillator phases based on user configuration
        if self.target_freq_sweep == "Secular":
            self.freq_update_arr = np.array([-1., 1., 0., 0., 0.])
        elif self.target_freq_sweep == "Probe":
            self.freq_update_arr = np.array([1., 1., 0., 0., 0.])
        if self.target_freq_sweep == "osc0":
            self.freq_update_arr = np.array([1., 0., 0., 0., 0.])
        elif self.target_freq_sweep == "osc1":
            self.freq_update_arr = np.array([0., 1., 0., 0., 0.])
        elif self.target_freq_sweep == "osc2":
            self.freq_update_arr = np.array([0., 0., 1., 0., 0.])
        elif self.target_freq_sweep == "osc3":
            self.freq_update_arr = np.array([0., 0., 0., 1., 0.])
        elif self.target_freq_sweep == "osc4":
            self.freq_update_arr = np.array([0., 0., 0., 0., 1.])

        # create pre-declared phase_offsets list as workaround for artiq stack memory issue
        # see: https://github.com/m-labs/artiq/issues/1520
        self.phase_offsets = np.zeros((len(self.phase_sweep_turns_list), 5))
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
            len(self.phase_heating_global_ch1_turns_list) *
            len(time_readout_mu_list),
        6), dtype=float)
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_experiment_list[:, [1, 2, -3, -2, -1, 0]] = np.stack(np.meshgrid(
            freq_carrier_hz_list,
            self.freq_sweep_hz_list,
            self.waveform_index_to_phase_sweep_turns,
            self.phase_heating_global_ch1_turns_list,
            time_readout_mu_list,
            freq_sideband_readout_ftw_list),
        -1).reshape(-1, 6)

        # randomize_config always enabled lol
        if self.randomize_config is True:
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
            if any((len(ampl_arg) != 5 for ampl_arg in ampl_arg_tuple)):
                raise ValueError("Error: phaser oscillator amplitude arrays must have length 5.")
            elif any((np.sum(ampl_arg) >= 100. for ampl_arg in ampl_arg_tuple)):
                raise ValueError("Error: phaser oscillator amplitudes must sum < 100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude arrays must be lists.")

        if isinstance(self.phase_osc_turns_list, list):
            if len(self.phase_osc_turns_list) != 5:
                raise ValueError("Error: phaser oscillator phase array must have length 5.")
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        if not isinstance(self.freq_osc_khz_list, list):
            raise ValueError("Error: phaser oscillator frequency array must be a list.")
        elif len(self.freq_osc_khz_list) != 5:
            raise ValueError("Error: phaser oscillator frequency array must have length 4.")

        if not isinstance(self.osc_ch1_sweep, list):
            raise ValueError("Error: phaser oscillator on/off must be a list.")
        elif len(self.osc_ch1_sweep) != 5:
            raise ValueError("Error: only 5 oscillators to change the phase of.")
        elif np.max(self.osc_ch1_sweep) > 1:
            raise ValueError("Error: can only turn oscillator on (1)/off (0)")
        elif np.max(self.osc_ch1_sweep) < 0:
            raise ValueError("Error: can only turn oscillator on (1)/off (0)")

        max_osc_freq_hz = (
                np.max(list(self.freq_sweep_khz_list)) * kHz +
                max(self.freq_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        min_osc_freq_hz = (
                np.min(list(self.freq_sweep_khz_list)) * kHz +
                max(self.freq_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_heating_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

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
        _osc_vals_ch0 = np.zeros((num_blocks, 5, 2), dtype=float)
        _osc_vals_ch1 = np.zeros((num_blocks, 5, 2), dtype=float)

        # update CH0 and CH1 sequences w/relevant amplitudes
        _osc_vals_ch0[:, :, 0] = np.array([self.ampl_ch0_stage_0, self.ampl_ch0_stage_1])
        _osc_vals_ch1[:, :, 0] = np.array([self.ampl_ch1_stage_0, self.ampl_ch1_stage_1])

        # phase track oscillator updates to account for 40ns sample period
        t_update_delay_s_list = np.array([0, 40e-9, 80e-9, 80e-9, 120e-9])
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_update_arr * np.mean(self.freq_sweep_hz_list)) *
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
                7)


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
        delay_mu(500000)

        # used to check_termination more frequently
        _loop_iter = 0

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep experiment configurations
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                carrier_freq_hz =   config_vals[1]
                freq_sweep_hz =     config_vals[2]
                phase_sweep_idx =   np.int32(config_vals[3])
                phase_ch1_turns =   config_vals[4]
                time_readout_mu =   np.int64(config_vals[5])

                # get corresponding phase and waveform ID from the index
                phase_sweep_turns = self.phase_sweep_turns_list[phase_sweep_idx]
                waveform_id = self.waveform_index_to_pulseshaper_id[phase_sweep_idx]
                self.core.break_realtime()

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_osc_base_hz_list + (freq_sweep_hz * self.freq_update_arr)
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    carrier_freq_hz - self.phaser_eggs.freq_center_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], freq_update_list[4]],
                    phase_ch1_turns
                )
                self.core.break_realtime()

                # set qubit readout frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_sb_readout, phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(50000)

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
                    phase_ch1_turns,
                    time_readout_mu
                )
                self.core.break_realtime()

                '''LOOP CLEANUP'''
                # resuscitate ion & detect deaths
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(counts)
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
        # should be backward-compatible with experiment which had no sub-repetitions
        try:
            sub_reps = self.sub_repetitions
        except Exception as e:
            sub_reps = 1

        # create lists for a ch1 sweep
        ch1_turns_sweep_list = []
        phonons_ch1_sweep_list = []

        # determine if a ch1 sweep occurred
        ch1_sweep_bool = len(self.phase_heating_global_ch1_turns_list) > 1
        phase_sweep_bool = len(self.phase_sweep_turns_list) > 1
        turns_sweep_bool = ch1_sweep_bool or phase_sweep_bool

        # print results
        print("\tResults - Heating:")

        # sweep over ch1_turns
        for ch1_turns in self.phase_heating_global_ch1_turns_list:

            # handle errors from data processing
            try:
                # convert dataset to array
                dataset_tmp = np.array(self.results)
                dataset = np.reshape(dataset_tmp[np.where(dataset_tmp[:, 5] == ch1_turns), :],
                                     (-1, dataset_tmp.shape[1]))

                ## determine scan type
                col_unique_vals = np.array([len(set(col)) for col in dataset.transpose()])
                # convert unique count to dataset index and order in decreasing value (excluding PMT counts)
                if np.argsort(-col_unique_vals)[0] != 1:
                    sorting_col_num = np.argsort(-col_unique_vals)[0]
                else:
                    sorting_col_num = np.argsort(-col_unique_vals)[1]

                # ensure we actually have a scan, and not some rubbish
                if col_unique_vals[sorting_col_num] <= 1:
                    continue

                ## convert results to sideband ratio
                ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
                                                                                               1, 0,
                                                                                               self.repetitions,
                                                                                               sub_reps)

                # get the phonons and instantiate the fitter class
                phonons = convert_ratios_to_coherent_phonons(ratios)
                fitter = fitSincGeneric()

                # format arguments for applet plotting
                if ch1_sweep_bool:
                    ccb_command = f'$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results'
                    group = 'plotting.results'
                    dataset_name = f'temp.plotting.results'
                    applet_name = f"CH 1 Ramsey"
                else:
                    ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results'
                    group = 'plotting.results'
                    dataset_name = f'temp.plotting.results'
                    applet_name = f"CH 1 Ramsey"

                ## process secular or carrier frequency sweep
                if sorting_col_num == 3 or sorting_col_num == 2:

                    fit_x = np.linspace(np.min(scanning_freq_MHz), np.max(scanning_freq_MHz),
                                        10 * len(scanning_freq_MHz))

                    # attempt to fit sinc to data
                    try:
                        # fit swept frequency (secular or carrier) and sidebands
                        fit_params_sweep, fit_err_sweep, _ = fitter.fit(scanning_freq_MHz, phonons)
                        fit_params_rsb, fit_err_rsb, _ = fitter.fit(scanning_freq_MHz, ave_rsb)
                        fit_params_bsb, fit_err_bsb, _ = fitter.fit(scanning_freq_MHz, ave_bsb, -1)

                        # get phonon from fit
                        phonon_n = fit_params_sweep[0]
                        # todo: implement phonon err
                        phonon_err = 0

                        # create arrays for plotting fits
                        fit_y_phonons = fitter.fit_func(fit_x, *fit_params_sweep)
                        fit_y_rsb = fitter.fit_func(fit_x, *fit_params_rsb)
                        fit_y_bsb = fitter.fit_func(fit_x, *fit_params_bsb)
                        # print results to log
                        print("\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_sweep[1] * 1e3,
                                                                          fit_err_sweep[1] * 1e3))

                        # save results to dataset manager for dynamic experiments
                        res_dj = [[phonon_n, phonon_err], [fit_params_sweep, fit_err_sweep]]
                        ch1_turns_sweep_list.append(ch1_turns)
                        phonons_ch1_sweep_list.append(phonon_n)

                        # save results to hdf5 as a dataset
                        if sorting_col_num == 3:
                            self.set_dataset('fit_params_secular', fit_params_sweep)
                            self.set_dataset('fit_err_secular', fit_err_sweep)
                        else:
                            self.set_dataset('fit_params_carrier', fit_params_sweep)
                            self.set_dataset('fit_err_carrier', fit_err_sweep)

                    # if fit fails then ignore plotting of fit by creating array of Nones
                    except Exception as e:
                        fit_y_phonons = [None] * len(fit_x)
                        fit_y_rsb = [None] * len(fit_x)
                        fit_y_bsb = [None] * len(fit_x)
                        res_dj = None

                    # format dictionary of results for plotting with applet
                    ccb_command += ' --num-subplots 3'
                    plotting_results = {'x': [scanning_freq_MHz, scanning_freq_MHz, scanning_freq_MHz],
                                        'y': [ave_rsb, ave_bsb, phonons],
                                        'errors': [std_rsb, std_bsb, [None] * len(std_rsb)],
                                        'fit_x': [fit_x, fit_x, fit_x],
                                        'fit_y': [fit_y_rsb, fit_y_bsb, fit_y_phonons],
                                        'subplot_x_labels': np.array(
                                            ['Frequency (MHz)', 'Frequency (MHz)', 'Frequency (MHz)']),
                                        'subplot_y_labels': np.array(
                                            ['D State Population', 'D State Population', 'Phonons']),
                                        'rid': self.scheduler.rid,
                                        }

                # process sideband readout sweep
                elif sorting_col_num == 0:
                    # todo: implement
                    phonon_err = 0
                    # get the sideband frequencies and prepare for plotting the fit
                    rsb_freqs_MHz, bsb_freqs_MHz, _ = extract_sidebands_freqs(scanning_freq_MHz)
                    fit_x_rsb = np.linspace(np.min(rsb_freqs_MHz), np.max(rsb_freqs_MHz), 1000)
                    fit_x_bsb = np.linspace(np.min(bsb_freqs_MHz), np.max(bsb_freqs_MHz), 1000)

                    try:
                        # try to fit the sidebands
                        fit_params_rsb, fit_err_rsb, fit_rsb = fitter.fit(rsb_freqs_MHz, ave_rsb)
                        fit_params_bsb, fit_err_bsb, fit_bsb = fitter.fit(bsb_freqs_MHz, ave_bsb)
                        fit_y_rsb = fitter.fit_func(fit_x_rsb, *fit_params_rsb)
                        fit_y_bsb = fitter.fit_func(fit_x_bsb, *fit_params_bsb)

                        # get the phonon number and update the list for graphing ch1 sweeps
                        phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
                        ch1_turns_sweep_list.append(ch1_turns)
                        phonons_ch1_sweep_list.append(phonon_n)

                        # save results to hdf5 as a dataset
                        self.set_dataset('fit_params_rsb', fit_params_rsb)
                        self.set_dataset('fit_params_bsb', fit_params_bsb)
                        self.set_dataset('fit_err_rsb', fit_err_rsb)
                        self.set_dataset('fit_err_bsb', fit_err_bsb)

                        # save results to dataset manager for dynamic experiments
                        res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]

                        # print results to log
                        print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2.,
                                                                  float(fit_err_rsb[1]) / 2.))
                        print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2.,
                                                                  float(fit_err_bsb[1]) / 2.))

                    except Exception as e:
                        print("Could not find fit the sidebands")
                        fit_y_rsb = [None] * len(fit_x_rsb)
                        fit_y_bsb = [None] * len(fit_x_bsb)
                        res_dj = None

                    # format dictionary for applet plotting
                    ccb_command += ' --num-subplots 2'
                    plotting_results = {'x': [rsb_freqs_MHz / 2, bsb_freqs_MHz / 2],
                                        'y': [ave_rsb, ave_bsb],
                                        'errors': [std_rsb, std_bsb],
                                        'fit_x': [fit_x_rsb / 2, fit_x_bsb / 2],
                                        'fit_y': [fit_y_rsb, fit_y_bsb],
                                        'subplot_x_labels': np.array(['AOM Frequency (MHz)', 'AOM Frequency (MHz)']),
                                        'subplot_y_labels': np.array(['D State Population', 'D State Population']),
                                        'rid': self.scheduler.rid,
                                        }

                else:
                    # if unknown quantity was scanned ignore analysis
                    plotting_results = {}
                    res_dj = None
                    group = None

                # record values in dataset
                self.set_dataset('temp.heating.results', res_dj, broadcast=True, persist=False,
                                 archive=False)
                self.set_dataset(dataset_name, pyon.encode(plotting_results), broadcast=True)
                self.set_dataset('temp.heating.rid', self.scheduler.rid, broadcast=True, persist=False,
                                 archive=False)

                # create applet
                self.ccb.issue("create_applet", applet_name, ccb_command, group=group)

            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))

        if ch1_sweep_bool:
            try:
                plotting_results = {'x': np.array(ch1_turns_sweep_list)[np.argsort(ch1_turns_sweep_list)],
                                    'y': np.array(phonons_ch1_sweep_list)[np.argsort(ch1_turns_sweep_list)],
                                    'subplot_x_labels': "Channel 1 Phase",
                                    'subplot_y_labels': 'Phonons',
                                    'rid': self.scheduler.rid
                                    }

                self.set_dataset('temp.plotting.results',
                                 pyon.encode(plotting_results), broadcast=True)

                self.ccb.issue("create_applet", f"CH1 Ramsey",
                               '$python -m LAX_exp.applets.plot_matplotlib '
                               'temp.plotting.results'
                               ' --num-subplots 1', group=['plotting', 'results'])


            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))
