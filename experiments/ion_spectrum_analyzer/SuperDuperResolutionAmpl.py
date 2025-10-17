from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from numpy import copy as np_copy
from numpy import array, int32, int64, zeros, arange, mean

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM,
    SidebandReadout, FockOverlap
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper, PULSESHAPER_MAX_WAVEFORMS


class SuperDuperResolutionAmpl(LAXExperiment, Experiment):
    """
    Experiment: Super Duper Resolution Ampl

    Investigate the QVSA (motional Raman) effect using phaser.
    Supports lots of easily configurable parameter scanning for phaser.
    Experiment name inspired by Sam Crary.
    """
    name = 'Super Duper Resolution Ampl'
    kernel_invariants = {
        # hardware values
        'att_phaser_mu', 'freq_global_offset_hz', 'freq_osc_base_hz_list',
        'waveform_index_to_compiled_wav', '_waveform_param_list',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence',
        'readout_subsequence', 'rescue_subsequence', 'fock_subsequence',
        'enable_RAP', 'enable_SBR',
        'spinecho_wizard', 'pulse_shaper',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',
        '_num_phaser_oscs', 'freq_readout_ftw_list',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=70, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config", BooleanValue(default=True),
                              tooltip="Completely randomizes the experiment sweep order.\n"
                                      "Needed b/c SDR-type exps are now used for EGGS searches.")
        self.setattr_argument("sub_repetitions", NumberValue(default=1, precision=0, step=1, min=1, max=500),
                              tooltip="Loop over the same config for a given number of sub_repetitions. "
                                      "Readout values are swept over for each sub_repetition.")
        self.setattr_argument("readout_type",   EnumerationValue(["SBR", "RAP", "RAP + SBR"], default="RAP"),
                              tooltip="SBR (Sideband Ratio): compares RSB and BSB amplitudes.\n"
                                      "RAP (Rapid Adiabatic Passage): Does RAP to measure fock state overlap.\n"
                                      "RAP + SBR: RAPs then SBRs (useful for e.g. fock overlap calibrations.")

        # todo: make this 5 and see if it's OK
        self._num_phaser_oscs = 4   # number of phaser oscillators in use

        # allocate relevant beam profiles
        self.profile_729_sb_readout =   0
        self.profile_729_SBC =          1
        self.profile_729_RAP =          2

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_sb_readout)
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)
        # extra argument for sideband readout (to enable rabiflop-type readout)
        self.setattr_argument("time_readout_us_list", Scannable(
                                                        default=[
                                                            ExplicitScan([26.11]),
                                                            RangeScan(0, 800, 200, randomize=True),
                                                        ],
                                                        global_min=1, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ),
                              group='sideband_readout',
                              tooltip="Sideband readout pulse times. Useful for rabiflop-type readout.")

        # fock state amplification
        self.fock_subsequence = FockOverlap(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=250,
            num_samples=200, pulse_shape="blackman"
        )

        # build experiment arguments
        self._build_arguments_freq_phase()
        self._build_arguments_waveform()
        self._build_arguments_PSK()

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def _build_arguments_freq_phase(self):
        """
        Set specific arguments for the phaser's frequency & phase.
        """
        _argstr = "SDR"  # create short string for argument grouping

        # phaser - frequency configuration
        self.setattr_argument("freq_phaser_carrier_mhz_list", Scannable(
                                                                        default=[
                                                                            ExplicitScan([86.]),
                                                                            CenterScan(86., 0.002, 0.0005, randomize=True),
                                                                        ],
                                                                        global_min=0.005, global_max=4800, global_step=1,
                                                                        unit="MHz", scale=1, precision=6
                                                                    ),
                              group="{}.freq".format(_argstr),
                              tooltip="Phaser output center frequency.\n"
                                      "Note: actual center frequency depends on the devices.phaser.freq_center_mhz dataset argument, "
                                      "which should be manually entered into the dataset manager by the user after "
                                      "configuring the TRF and NCO via e.g. the phaser_configure tool.\n"
                                      "Ensure all values are set correctly.")
        self.setattr_argument("freq_sweep_arr", PYONValue([-1., 1., 0., 0.]),
                              group="{}.freq".format(_argstr),
                              tooltip="Defines how oscillator freqs should be scaled for values in freq_osc_sweep_khz_list.\n"
                                      "Indices of freq_sweep_arr correspond to the oscillator number. "
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the freq value, and osc_1 by -1x the freq value, with the rest untouched.\n"
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("freq_osc_sweep_khz_list",    Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        CenterScan(0., 4, 0.1, randomize=True),
                                                                        RangeScan(0., 100.0, 26, randomize=True),
                                                                    ],
                                                                    global_min=-10000, global_max=10000, global_step=10,
                                                                    unit="kHz", scale=1, precision=6
                                                                ),
                              group="{}.freq".format(_argstr),
                              tooltip="Frequency sweep applied via the phaser oscillators.\n"
                                      "Values for each oscillator are adjusted by the array in freq_sweep_arr.")

        # phaser - phase configuration
        self.setattr_argument("phase_sweep_arr", PYONValue([1., 0., 0., 0.]),
                              group="{}.phase".format(_argstr),
                              tooltip="Defines how oscillator phases should be adjusted for each value in phase_osc_sweep_turns_list. "
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the phase value, and osc_1 by -1x the phase value, with the rest untouched. "
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("phase_osc_sweep_turns_list", Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        RangeScan(0, 1.0, 26, randomize=True),
                                                                    ],
                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                    unit="turns", scale=1, precision=5
                                                                ),
                              group="{}.phase".format(_argstr),
                              tooltip="Phase sweep applied via the phaser oscillators.\n"
                                      "Values for each oscillator are adjusted by the array in phase_sweep_arr.")
        self.setattr_argument("phase_global_ch1_turns_list",  Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        RangeScan(0, 1.0, 26, randomize=True),
                                                                    ],
                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                    unit="turns", scale=1, precision=5
                                                                ),
                              group="{}.phase".format(_argstr),
                              tooltip="Sets a global CH1 phase via the DUC.\n"
                                      "Note: the devices.phaser.phas_ch1_inherent_turns dataset argument is overridden "
                                      "in this experiment.")

    def _build_arguments_waveform(self):
        """
        Set specific arguments for the phaser waveform.
        """
        _argstr = "SDR"  # create short string for argument grouping

        # phaser - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False),
                              group='{}.shape'.format(_argstr),
                              tooltip="Applies pulse shaping to the edges of the phaser pulse.\n"
                                      "Note: pulse shaping is applied to each constituent PSK block.")
        self.setattr_argument("type_pulse_shape",       EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.shape'.format(_argstr),
                              tooltip="Pulse shape type to be used.")
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Time constant of the pulse shape. This is used for both the pulse rollon AND rolloff.\n"
                                      "e.g. a 1ms main pulse time with 100us time_pulse_shape_rolloff_us will result in a 1ms + 2*100us = 1.2ms total pulse time.\n"
                                      "All constituent PSK blocks will have this pulse time applied.\n"
                                      "Note: DMA issues limit the total number of samples (i.e. time_pulse_shape_rolloff_us * freq_pulse_shape_sample_khz).")
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Sample rate used for pulse shaping.\n"
                                      "This value is inexact and is fixed at multiples of the phaser oscillator update "
                                      "rate (i.e. 40ns) times the number of oscillators in use.")

        # phaser - waveform - general
        self.setattr_argument("time_heating_us",   NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.wav".format(_argstr),
                              tooltip="Total MAIN pulse time per phaser pulse. "
                                      "This time is split among all PSK blocks and does not include any PSK delays.\n"
                                      "IMPORTANT NOTE: pulse shaping times are IN ADDITION to time_heating us, and ALL PSK blocks are pulse shaped.\n"
                                      "e.g. 1ms time_heating_us with 21 PSKs and 100us time_pulse_shape_rolloff_us results in an actual pulse time of "
                                      "1ms + 21 * (100us * 2) = 5.2ms.")
        self.setattr_argument("att_phaser_db",    NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.wav".format(_argstr),
                              tooltip="Phaser attenuation to be used for both CH0 and CH1.")
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.wav".format(_argstr),
                              tooltip="Apply a frequency offset via the phaser oscillators to avoid any DUC/NCO/TRF output spurs.\n"
                                      "Range is limited by the phaser oscillator freq range, i.e. [-10, 10] MHz (includes the frequencies in freq_osc_khz_list).")
        self.setattr_argument("freq_osc_khz_list",      PYONValue([-702.687, 702.687, 0., 0.]),
                              group="{}.wav".format(_argstr),
                              tooltip="Phaser oscillator frequencies.")
        self.setattr_argument("ampl_osc_frac_list",     PYONValue([40., 40., 10., 0.]),
                              group="{}.wav".format(_argstr),
                              tooltip="Phaser oscillator amplitudes. Applied to both CH0 and CH1.\n"
                                      "Note: CH1 amplitudes will be scaled by the amplitude scaling factors in devices.phaser.ch1.ampl_ch1_osc_scale_arr.")
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0.]),
                              group="{}.wav".format(_argstr),
                              tooltip="Relative phases between each phaser oscillator. Applied on both CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0.5, 0.5, 0.5]),
                              group="{}.wav".format(_argstr),
                              tooltip="Sets the relative CH1 phase via the phaser oscillators.")


        # phase - waveform - sweep
        self.setattr_argument("enable_wav_scale",   BooleanValue(default=False),
                              group='{}.wav_sweep'.format(_argstr),
                              tooltip="Scale relevant waveform parameters (e.g. amplitude, time).")
        self.setattr_argument("target_wav_scale",   EnumerationValue(['Amplitude', 'Time (Total)', 'Time (Main)', 'Time (Shape)'], default='Amplitude'),
                              group='{}.wav_sweep'.format(_argstr),
                              tooltip="Select the waveform parameter to be scaled.\n"
                                      "Amplitude: scales the amplitudes of all oscillators.\n"
                                      "Time (Total): scales the total time (i.e. main pulse AND pulse-shaped edges).\n"
                                      "Time (Main): scales only the main pulse time (i.e. excludes the pulse-shaped edges).\n"
                                      "Time (Shape): scales only the time constant for the pulse-shaped edges (i.e. excludes the main pulse).\n")
        self.setattr_argument("wav_osc_scale_list", Scannable(
                                                        default=[
                                                            ExplicitScan([1.]),
                                                            CenterScan(0.5, 0.5, 0.02, randomize=True),
                                                            RangeScan(0., 1.0, 26, randomize=True),
                                                        ],
                                                        global_min=0., global_max=1., global_step=0.05,
                                                        scale=1, precision=5
                                                    ),
                              group="{}.wav_sweep".format(_argstr),
                              tooltip="Waveform scaling factor applied to ampl_osc_frac_list.\n"
                                      "Enables the amplitude to be scanned when e.g. optimizing powers.\n"
                                      "Scaling factor is applied equally to BOTH phaser output channels.")

    def _build_arguments_PSK(self):
        """
        Set specific arguments for PSK scheduling.
        """
        _argstr = "SDR"  # create short string for argument grouping

        # phaser - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",  BooleanValue(default=False), group="{}.psk".format(_argstr),
                              tooltip="Enable PSK-ing: break the main pulse into individual blocks with different phases.\n"
                                      "Number of PSKs is determined by number of phases in phase_osc<x>_psk_turns. "
                                      "All oscillator PSK schedules must have same length.")
        self.setattr_argument("phase_osc0_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc0.")
        self.setattr_argument("phase_osc1_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc1.")
        self.setattr_argument("phase_osc2_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc2.")
        self.setattr_argument("phase_osc3_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc3.")
        self.setattr_argument("phase_osc4_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc4.")

        # phaser - waveform - PSK (i.e. ramsey-style) delay
        self.setattr_argument("enable_psk_delay", BooleanValue(default=False),
                              group="{}.psk".format(_argstr),
                              tooltip="Add a delay between PSK pulses where oscillator amplitudes are set to 0. "
                                      "Can be used to create e.g. a Ramsey or DD-type pulse sequence.\n"
                                      "Requires enable_phase_shift_keying to be enabled; otherwise, does nothing.\n"
                                      "Note: prepare/cleanup methods (e.g. set phaser atts, set ext switch) are not called for the delay.")
        self.setattr_argument("time_psk_delay_us_list", Scannable(
                                                      default=[
                                                          ExplicitScan([1000.]),
                                                          RangeScan(10, 1000, 20, randomize=True),
                                                      ],
                                                      global_min=1, global_max=100000, global_step=1,
                                                      unit="us", scale=1, precision=5
                                                  ),
                              group="{}.psk".format(_argstr),
                              tooltip="Delay time (in us) delay between PSK pulses. Used for ramsey-ing.\n"
                                      "Note: enable_phase_shift_keying AND enable_psk_delay must be True.")

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''
        GENERAL SETUP
        '''
        self._prepare_argument_checks()

        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_osc_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, array(self.phase_osc_ch1_offset_turns))

        # configure readout method
        self.enable_RAP, self.enable_SBR = (False, False) # default to false for both
        if 'RAP' in self.readout_type:
            self.enable_RAP = True
            self.freq_readout_ftw_list = array([-1], dtype=int32)
            time_readout_mu_list = [256]
        # note: SBR check has to be 2nd, since otherwise RAP overrides SBR's freq and time list
        if 'SBR' in self.readout_type:
            self.enable_SBR = True
            self.freq_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
            time_readout_mu_list = [self.core.seconds_to_mu(time_us * us) for time_us in self.time_readout_us_list]
        if not any(kw in self.readout_type for kw in ('RAP', 'SBR')):
            raise ValueError("Invalid readout type. Must be one of (SBR, RAP, RAP + SBR).")


        '''
        HARDWARE VALUES - CONFIG
        '''
        # convert values to convenience units
        self.att_phaser_mu = att_to_mu(self.att_phaser_db * dB)
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # format build arguments as numpy arrays with appropriate units
        freq_phaser_carrier_hz_list = array(list(self.freq_phaser_carrier_mhz_list)) * MHz
        freq_osc_sweep_hz_list = array(list(self.freq_osc_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_hz
        self.phase_osc_sweep_turns_list = list(self.phase_osc_sweep_turns_list)
        self.freq_sweep_arr = array(self.freq_sweep_arr, dtype=float)
        self.phase_sweep_arr = array(self.phase_sweep_arr, dtype=float)


        '''
        EXPERIMENT CONFIGURATION
        '''
        # collate waveform sweep parameters into a list for later batch compilation
        # todo: ensure delay times only in 320ns steps
        if self.enable_phase_shift_keying and self.enable_psk_delay:
            time_psk_delay_us_list_dj = list(self.time_psk_delay_us_list)
        else: time_psk_delay_us_list_dj = [0]

        if self.enable_wav_scale: wav_osc_scale_list = list(self.wav_osc_scale_list)
        else: wav_osc_scale_list = [1.]

        # create waveform parameter sweep config (for use by _prepare_waveform)
        self._waveform_param_list = create_experiment_config(
            self.phase_osc_sweep_turns_list,
            time_psk_delay_us_list_dj,
            wav_osc_scale_list,
            shuffle_config=False, config_type=float
        )
        waveform_num_list = arange(len(self._waveform_param_list))

        # create actual experiment config
        self.config_experiment_list = create_experiment_config(
            freq_phaser_carrier_hz_list, freq_osc_sweep_hz_list,
            waveform_num_list, list(self.phase_global_ch1_turns_list),
            time_readout_mu_list,
            config_type=float, shuffle_config=self.randomize_config
        )

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''
        CHECK PHASER BASE OSC CONFIG
        '''
        # check that phaser oscillator amplitude config is valid
        if ((not isinstance(self.ampl_osc_frac_list, list)) or
                (len(self.ampl_osc_frac_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator amplitude array must be list of length {:d}.".format(self._num_phaser_oscs))
        elif sum(self.ampl_osc_frac_list) >= 100.:
            raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
        # todo: account for wav_osc_scale_list

        # check that phaser oscillator phase arrays are valid
        if ((not isinstance(self.phase_osc_turns_list, list)) or
                (len(self.phase_osc_turns_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator phase array must be list of length {:d}.".format(self._num_phaser_oscs))

        # check that phaser oscillator frequencies are valid
        if ((not isinstance(self.freq_osc_khz_list, list)) or
                (len(self.freq_osc_khz_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator frequency array must be list of length {:d}.".format(self._num_phaser_oscs))
        max_osc_freq_hz = (max(list(self.freq_osc_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        min_osc_freq_hz = (min(list(self.freq_osc_sweep_khz_list)) * kHz +
                           min(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = array(list(self.freq_phaser_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 300. * MHz) or (phaser_carrier_lower_dev_hz >= 300. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")


        '''
        CHECK PHASER WAVEFORM CONFIG
        '''
        # check that PSK schedule is valid
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (not isinstance(psk_schedule, list)) or (len(psk_schedule) != num_psk_blocks)
            for psk_schedule in (
                self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                self.phase_osc2_psk_turns, self.phase_osc3_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule: all PSK schedules must be of same length.")

        # ensure that sweep targets are lists of appropriate length
        if not (isinstance(self.phase_sweep_arr, list) and (len(self.phase_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid phase_sweep_arr: {:}.\nphase_sweep_arr must be list of length {:d}.".format(
                self.phase_sweep_arr, self._num_phaser_oscs))
        if not (isinstance(self.freq_sweep_arr, list) and (len(self.freq_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid freq_sweep_arr: {:}.\nfreq_sweep_arr must be list of length {:d}.".format(
                self.freq_sweep_arr, self._num_phaser_oscs))

        # check that waveforms are not too many/not sweeping too hard
        num_waveforms_to_record = (len(list(self.phase_osc_sweep_turns_list)) *
                                   len(list(self.time_psk_delay_us_list)) *
                                   len(list(self.wav_osc_scale_list)))
        if num_waveforms_to_record > PULSESHAPER_MAX_WAVEFORMS:
            raise ValueError("Too many waveforms to record ({:d}) - must be fewer than {:d}.\n"
                             "Reduce length of any of [phase_osc_sweep_turns_list, "
                             "time_psk_delay_us_list, wav_osc_scale_list].".format(
                num_waveforms_to_record, PULSESHAPER_MAX_WAVEFORMS))

        # todo: check scaling aligns with pulse shaping

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA/phaser pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''
        PREPARE WAVEFORM COMPILATION
        '''
        # instantiate relevant global variables
        self.waveform_index_to_compiled_wav = list() # store compiled waveform values
        # note: waveform_index_to_pulseshaper_id is NOT kernel_invariant b/c gets updated in phaser_record
        self.waveform_index_to_pulseshaper_id = zeros(len(self._waveform_param_list), dtype=int32) # store waveform ID linked to DMA sequence

        # calculate block timings and scale amplitudes for ramsey-ing
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        if self.enable_phase_shift_keying:

            # case: spin-echo style (PSK w/ramsey delay)
            if self.enable_psk_delay:
                # note: don't create block_time_list_us here b/c we need to scan ramsey delay time
                block_ampl_scale_list = riffle([1] * num_psk_blocks, [0] * (num_psk_blocks - 1))

            # case: PSK-style (PSK only, no ramsey delay)
            else:
                block_time_list_us = [self.time_heating_us / num_psk_blocks] * num_psk_blocks
                block_ampl_scale_list = [1] * num_psk_blocks

        # case: rabi-style (no PSK, no ramsey delay)
        else:
            block_time_list_us = [self.time_heating_us]
            block_ampl_scale_list = [1]


        '''
        DESIGN WAVEFORM SEQUENCE
        '''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = zeros((len(block_time_list_us), self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = array(self.ampl_osc_frac_list)  # set oscillator amplitudes
        # use block_ampl_scale_list to disable output (by setting amplitudes to zero) during delay blocks
        _osc_vals_blocks[:, :, 0] *= array([block_ampl_scale_list]).transpose()

        # set oscillator phases and account for oscillator update delays
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        freq_osc_sweep_avg_hz = mean(list(self.freq_osc_sweep_khz_list)) * kHz
        t_update_delay_s_list = array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_sweep_arr * freq_osc_sweep_avg_hz) *
                t_update_delay_s_list
        )
        _osc_vals_blocks[:, :, 1] += array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # if we have delays between PSK blocks, use ::2 since we want to only update non-delay blocks (and not delay blocks)
                _osc_vals_blocks[::2, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                      self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                      self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()
            else:
                # otherwise, apply PSK schedule to all blocks
                _osc_vals_blocks[:, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                    self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                    self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()


        '''
        CREATE PARAMETER-SPECIFIC WAVEFORMS
        '''
        # record phaser waveforms - one for each phase
        for waveform_params in self._waveform_param_list:
            # extract waveform params
            phase_sweep_turns = waveform_params[0]
            time_us_delay = waveform_params[1]
            wav_osc_scale_val = waveform_params[2]


            '''
            PREPARE LOCAL PARAMETER-SPECIFIC CONFIGS 
            '''
            # create local working copy of block_time_list_us
            # case: spin-echo style (PSK w/ramsey delay) - create block_time_list_us here
            #   since the target delay time is a scannable we can only set here
            if self.enable_phase_shift_keying and self.enable_psk_delay:
                _block_time_list_us_local = riffle([self.time_heating_us / num_psk_blocks] * num_psk_blocks,
                                                   [time_us_delay] * (num_psk_blocks - 1))
            # case: all others - simply copy block_time_list_us
            else: _block_time_list_us_local = np_copy(block_time_list_us)

            # create local working copy of _osc_vals_blocks and update with waveform parameters
            # note: no need to deep copy b/c it's filled w/immutables
            _osc_vals_blocks_local = np_copy(_osc_vals_blocks)
            _osc_vals_blocks_local[:, :, 1] += self.phase_sweep_arr * phase_sweep_turns # apply phase sweep


            # waveform sweep: scale waveform parameter according to target_wav_scale
            _time_pulse_shape_rolloff_local = self.time_pulse_shape_rolloff_us
            if self.enable_wav_scale:

                # ensure valid wav_osc_scale_val for scaling
                if 0 <= wav_osc_scale_val <= 1.:

                    # case - "Amplitude": scale all oscillator amplitues
                    if self.target_wav_scale == "Amplitude":
                        _osc_vals_blocks_local[:, :, 0] *= wav_osc_scale_val

                    # case - ("Time (Total)", "Time (Shape)"): scale pulse-shaped edges
                    if any(kw == self.target_wav_scale for kw in ('Time (Total)', 'Time (Shape)')):
                        _time_pulse_shape_rolloff_local *= wav_osc_scale_val    # scale pulse-shaped rolloff times

                    # case - ("Time (Total)", "Time (Main)"): scale main pulse time
                    if any(kw == self.target_wav_scale for kw in ('Time (Total)', 'Time (Main)')):

                        # case - spin-echo style: scale only non-delay blocks
                        if self.enable_phase_shift_keying and self.enable_psk_delay:
                            _block_time_list_us_local = [
                                time_us * wav_osc_scale_val
                                if i % 2 == 0 else time_us
                                for i, time_us in enumerate(_block_time_list_us_local)
                            ]

                        # case - all other styles: scale all blocks
                        else:
                            _block_time_list_us_local = [
                                time_us * wav_osc_scale_val
                                for time_us in _block_time_list_us_local
                            ]

                else:
                    raise ValueError("Invalid wav_osc_scale_val - must be in [0., 1.].")


            '''
            CREATE WAVEFORM SEQUENCE DICT & COMPILE 
            '''
            # specify sequence as a list of blocks, where each block is a dict
            # note: have to instantiate _sequence_blocks_local locally each time b/c dicts aren't deep copied
            _sequence_blocks_local = [
                {
                    "oscillator_parameters": _osc_vals_blocks_local[_idx_block],
                    "config": {
                        "time_us": _block_time_list_us_local[_idx_block],
                        # don't pulse shape for delay blocks lmao
                        "pulse_shaping": self.enable_pulse_shaping and (block_ampl_scale_list[_idx_block] != 0),
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping,
                            "pulse_shape_falling": self.enable_pulse_shaping,
                            "sample_rate_khz": self.freq_pulse_shape_sample_khz,
                            "rolloff_time_us": _time_pulse_shape_rolloff_local
                        }
                    }
                } for _idx_block in range(len(_block_time_list_us_local))
            ]
            # compile waveform to numerical values and store in holder
            self.waveform_index_to_compiled_wav.append(
                self.spinecho_wizard.compile_waveform(_sequence_blocks_local))

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions *
                len(self.config_experiment_list) * len(self.freq_readout_ftw_list),
                9)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # record fock subsequences
        with self.core_dma.record('_FOCK_GEN_HANDLE'):
            self.fock_subsequence.run_fock_generate()
        with self.core_dma.record('_FOCK_READ_HANDLE'):
            self.fock_subsequence.run_fock_read()

        # record phaser waveforms
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage during configuration
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # load extra DMA handles
        self.pulse_shaper.waveform_load()
        fock_gen_handle = self.core_dma.get_handle('_FOCK_GEN_HANDLE')
        fock_read_handle = self.core_dma.get_handle('_FOCK_READ_HANDLE')

        _loop_iter = 0  # used to check_termination more frequently

        '''MAIN LOOP'''
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''
                CONFIGURE
                '''
                # extract values from config list
                carrier_freq_hz =   config_vals[0]
                freq_sweep_hz =     config_vals[1]
                waveform_num =      int32(config_vals[2])
                phase_ch1_turns =   config_vals[3]
                time_readout_mu =   int64(config_vals[4])

                # get corresponding waveform parameters and pulseshaper ID from the index
                waveform_params = self._waveform_param_list[waveform_num]
                phaser_waveform = self.waveform_index_to_pulseshaper_id[waveform_num]

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_osc_base_hz_list + freq_sweep_hz * self.freq_sweep_arr

                self.core.break_realtime()
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    carrier_freq_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], 0.],
                    # CH1 phase (via DUC)
                    phase_ch1_turns
                )


                '''
                SUB-REPETITION IMPLEMENTATION
                '''
                for _subrep_num in range(self.sub_repetitions):
                    # sub-repetitions requires we run a bunch of sub_reps at the same config
                    # however, we still want to sweep the readout parameters
                    for freq_readout_ftw in self.freq_readout_ftw_list:

                        # configure readout
                        self.core.break_realtime()
                        if self.enable_SBR:
                            self.qubit.set_mu(freq_readout_ftw,
                                              asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                              profile=self.profile_729_sb_readout,
                                              phase_mode=PHASE_MODE_CONTINUOUS)
                            delay_mu(25000)

                        '''STATE PREPARATION'''
                        self.initialize_subsequence.run_dma()
                        self.sidebandcool_subsequence.run_dma()
                        self.core_dma.playback_handle(fock_gen_handle) # generate higher fock states

                        '''QVSA PULSE'''
                        self.phaser_run(phaser_waveform)

                        '''READOUT'''
                        if self.enable_RAP:
                            self.core_dma.playback_handle(fock_read_handle)
                        if self.enable_SBR:
                            self.sidebandreadout_subsequence.run_time(time_readout_mu)
                        self.readout_subsequence.run_dma()
                        self.rescue_subsequence.resuscitate()

                        '''LOOP CLEANUP'''
                        counts = self.readout_subsequence.fetch_count()
                        self.rescue_subsequence.detect_death(counts)
                        self.update_results(
                            freq_readout_ftw,
                            counts,
                            carrier_freq_hz,
                            freq_sweep_hz,
                            waveform_params[0], # note: manually expand waveform_params b/c no variadics in kernel
                            waveform_params[1],
                            waveform_params[2],
                            phase_ch1_turns,
                            time_readout_mu
                        )

                        # check termination more frequently in case reps are low
                        if _loop_iter % 50 == 0:
                            self.check_termination()
                        _loop_iter += 1

            # rescue ion as needed & check termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    HELPER FUNCTIONS - PHASER
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main phaser pulse together with supporting functionality.
        :param waveform_id: the ID of the waveform to run.
        """
        # QVSA - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_phaser_mu, self.att_phaser_mu)

        # QVSA - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # QVSA - STOP
        # stop all output & clean up hardware (e.g. phaser amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each waveform parameter
        for _wav_idx in range(len(self._waveform_param_list)):
            # get waveform for given parameters
            # note: use sync RPC to reduce significant overhead of direct data transfer
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self._get_compiled_waveform(_wav_idx)

            # record phaser pulse sequence and save returned waveform ID
            # note: no need to add slack b/c waveform_record does it for us
            self.waveform_index_to_pulseshaper_id[_wav_idx] = self.pulse_shaper.waveform_record(
                _wav_data_ampl, _wav_data_phas, _wav_data_time
            )

    @rpc
    def _get_compiled_waveform(self, wav_idx: TInt32) -> TTuple([TArray(TFloat, 2),
                                                                 TArray(TFloat, 2),
                                                                 TArray(TInt64, 1)]):
        """
        Return compiled waveform values.
        By returning the large waveform arrays via RPC, we avoid all-at-once large data transfers,
            speeding up experiment compilation and transfer to Kasli.
        :param wav_idx: the index of the compiled waveform to retrieve.
        :return: the compiled waveform corresponding to the waveform index.
        """
        return self.waveform_index_to_compiled_wav[wav_idx]

