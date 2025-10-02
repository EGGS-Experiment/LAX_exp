from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_TRACKING

from numpy import copy as np_copy
from numpy import array, int32, zeros, arange, mean

from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper, PULSESHAPER_MAX_WAVEFORMS
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes


class PSRSB_RDX(LAXExperiment, Experiment):
    """
    Experiment: Phase-Sensitive Red Sideband (PSRSB) RDX

    Use QVSA to generate a phase-sensitive displacement, then apply the PSRSB technique for readout.
    """
    name = 'PSRSB RDX'
    kernel_invariants = {
        # hardware values - phaser
        'att_phaser_mu', 'freq_osc_sweep_hz_list', 'freq_phaser_carrier_hz',
        'freq_osc_base_hz_list', 'waveform_index_to_compiled_wav', '_waveform_param_list',

        # hardware values - PSRSB
        'att_qubit_mu', 'ampl_psrsb_rsb_asf', 'ampl_psrsb_carr_asf', 'time_psrsb_rsb_mu', 'time_psrsb_carr_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'spinecho_wizard', 'pulse_shaper',

        # configs
        'profile_729_PSRSB', 'profile_729_SBC', 'config_experiment_list', '_num_phaser_oscs',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=70, precision=0, step=1, min=1, max=100000))
        self._num_phaser_oscs = 4   # number of phaser oscillators in use

        # allocate relevant beam profiles
        self.profile_729_PSRSB =    0
        self.profile_729_SBC =      1

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # set build arguments
        self._build_arguments_freq_phase()
        self._build_arguments_waveform()
        self._build_arguments_PSK()
        self._build_arguments_PSRSB()

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def _build_arguments_freq_phase(self):
        """
        Set specific arguments for the phaser's frequency & phase.
        """
        _argstr = "SDR"  # create short string for argument grouping

        # phaser - frequency configuration
        self.setattr_argument("freq_phaser_carrier_mhz",    NumberValue(default=86., precision=7, step=1, min=0.001, max=4800, unit="MHz", scale=1.),
                              group="{}.global".format(_argstr))
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
                                      "Note: the eggs.phas_ch1_inherent_turns dataset argument is overridden "
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
                              group="{}.waveform".format(_argstr),
                              tooltip="Total MAIN pulse time per phaser pulse. "
                                      "This time is split among all PSK blocks and does not include any PSK delays.\n"
                                      "IMPORTANT NOTE: pulse shaping times are IN ADDITION to time_heating us, and ALL PSK blocks are pulse shaped.\n"
                                      "e.g. 1ms time_heating_us with 21 PSKs and 100us time_pulse_shape_rolloff_us results in an actual pulse time of "
                                      "1ms + 21 * (100us * 2) = 5.2ms.")
        self.setattr_argument("att_phaser_db",    NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser attenuation to be used for both CH0 and CH1.")
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Apply a frequency offset via the phaser oscillators to avoid any DUC/NCO/TRF output spurs.\n"
                                      "Range is limited by the phaser oscillator freq range, i.e. [-10, 10] MHz (includes the frequencies in freq_osc_khz_list).")
        self.setattr_argument("freq_osc_khz_list",      PYONValue([-702.687, 702.687, 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator frequencies.")
        self.setattr_argument("ampl_osc_frac_list",     PYONValue([40., 40., 10., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator amplitudes. Applied to both CH0 and CH1.\n"
                                      "Note: CH1 amplitudes will be scaled by the amplitude scaling factors in eggs.ch1.ampl_ch1_osc_scale_arr.")
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Relative phases between each phaser oscillator. Applied on both CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0.5, 0.5, 0.5]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Sets the relative CH1 phase via the phaser oscillators.")

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
        self.setattr_argument("enable_psk_delay", BooleanValue(default=False), group="{}.psk".format(_argstr),
                              tooltip="Add a delay between PSK pulses where oscillator amplitudes are set to 0. "
                                      "Can be used to create e.g. a Ramsey or DD-type pulse sequence. "
                                      "Requires enable_phase_shift_keying to be enabled; otherwise, does nothing.\n"
                                      "Note: prepare/cleanup methods (e.g. set phaser atts, set ext switch) are not called for the delay.")
        self.setattr_argument("time_psk_delay_us", NumberValue(default=100, precision=3, step=100, min=0.2, max=1000000, unit="us", scale=1.),
                              group="{}.psk".format(_argstr),
                              tooltip="Delay time (in us) delay between PSK pulses. Used for ramsey-ing.\n"
                                      "Note: enable_phase_shift_keying AND enable_psk_delay must be True.")

    def _build_arguments_PSRSB(self):
        """
        Set specific arguments for PSRSB beam configuration.
        """
        _argstr = "PSRSB"    # string to use for arguments

        # PSRSB - general
        self.setattr_argument("att_qubit_db", NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group='{}'.format(_argstr))

        # PSRSB - RSB pulse
        self.setattr_argument("ampl_psrsb_rsb_pct", NumberValue(default=50, precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group='{}.RSB'.format(_argstr))
        self.setattr_argument("time_psrsb_rsb_us",  NumberValue(default=27.51, precision=3, step=5, min=1, max=10000000, scale=1., unit="us"),
                              group='{}.RSB'.format(_argstr))
        self.setattr_argument("freq_psrsb_rsb_mhz_list",    Scannable(
                                                                default=[
                                                                    CenterScan(100.7471, 0.01, 0.0001, randomize=True),
                                                                    ExplicitScan([100.7471]),
                                                                ],
                                                                global_min=60., global_max=200., global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group='{}.RSB'.format(_argstr))
        self.setattr_argument("phas_psrsb_rsb_turns_list",  Scannable(
                                                                default=[
                                                                    RangeScan(0, 1.0, 21, randomize=True),
                                                                    ExplicitScan([0.372]),
                                                                ],
                                                                global_min=-1.0, global_max=1.0, global_step=0.1,
                                                                unit="turns", scale=1, precision=3
                                                            ), group='{}.RSB'.format(_argstr))

        # PSRSB - carrier pulse
        self.setattr_argument("ampl_psrsb_carr_pct",         NumberValue(default=50, precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group='{}.carr'.format(_argstr))
        self.setattr_argument("time_psrsb_carr_us",          NumberValue(default=1.21, precision=3, step=1, min=1, max=10000000, scale=1., unit="us"),
                              group='{}.carr'.format(_argstr))
        self.setattr_argument("freq_psrsb_carr_mhz_list",    Scannable(
                                                                default=[
                                                                    ExplicitScan([101.0981]),
                                                                    CenterScan(101.0981, 0.01, 0.0001, randomize=True),
                                                                ],
                                                                global_min=60., global_max=200., global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group='{}.carr'.format(_argstr))
        self.setattr_argument("phas_psrsb_carr_turns_list",  Scannable(
                                                                default=[
                                                                    RangeScan(0, 1.0, 21, randomize=True),
                                                                    ExplicitScan([0.372]),
                                                                ],
                                                                global_min=-1.0, global_max=1.0, global_step=0.1,
                                                                unit="turns", scale=1, precision=3
                                                            ), group='{}.carr'.format(_argstr))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE & VALIDATE INPUTS'''
        self._prepare_argument_checks()


        '''HARDWARE VALUES - PHASER'''
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_osc_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, array(self.phase_osc_ch1_offset_turns))

        # convert values to convenience units
        self.att_phaser_mu = att_to_mu(self.att_phaser_db * dB)
        self.freq_phaser_carrier_hz = (self.freq_phaser_carrier_mhz - self.freq_global_offset_mhz) * MHz

        # format build arguments as numpy arrays with appropriate units
        self.freq_osc_sweep_hz_list = array(list(self.freq_osc_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_mhz * MHz
        self.phase_osc_sweep_turns_list = list(self.phase_osc_sweep_turns_list)
        self.freq_sweep_arr = array(self.freq_sweep_arr, dtype=float)
        self.phase_sweep_arr = array(self.phase_sweep_arr, dtype=float)


        '''HARDWARE VALUES - PSRSB'''
        # PSRSB DDS waveform params
        self.att_qubit_mu = self.qubit.cpld.att_to_mu(self.att_qubit_db * dB)

        self.ampl_psrsb_rsb_asf =   self.qubit.amplitude_to_asf(self.ampl_psrsb_rsb_pct / 100.)
        self.ampl_psrsb_carr_asf =  self.qubit.amplitude_to_asf(self.ampl_psrsb_carr_pct / 100.)

        self.time_psrsb_rsb_mu =        self.core.seconds_to_mu(self.time_psrsb_rsb_us * us)
        self.time_psrsb_carr_mu =    self.core.seconds_to_mu(self.time_psrsb_carr_us * us)

        # PSRSB waveform scans
        freq_psrsb_rsb_ftw_list =   [self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                     for freq_mhz in list(self.freq_psrsb_rsb_mhz_list)]
        freq_psrsb_carr_ftw_list =  [self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                     for freq_mhz in list(self.freq_psrsb_carr_mhz_list)]
        phas_psrsb_rsb_pow_list =   [self.qubit.turns_to_pow(phas_turns)
                                     for phas_turns in list(self.phas_psrsb_rsb_turns_list)]
        phas_psrsb_carr_pow_list =  [self.qubit.turns_to_pow(phas_turns)
                                     for phas_turns in list(self.phas_psrsb_carr_turns_list)]

        '''EXPERIMENT CONFIGURATION'''
        # collate waveform sweep parameters into a list for later batch compilation
        # note: we do things in this roundabout way b/c want to support large phaser waveform sweeps
        #   for further notes, compare to SDR exp
        self._waveform_param_list = create_experiment_config(self.phase_osc_sweep_turns_list,
                                                             shuffle_config=False,
                                                             config_type=float)
        waveform_num_list = arange(len(self._waveform_param_list))

        # create experiment config
        self.config_experiment_list = create_experiment_config(
            self.freq_osc_sweep_hz_list,
            waveform_num_list,
            list(self.phase_global_ch1_turns_list),
            freq_psrsb_rsb_ftw_list,
            freq_psrsb_carr_ftw_list,
            phas_psrsb_rsb_pow_list,
            phas_psrsb_carr_pow_list,
            config_type=float,
            shuffle_config=True
        )

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''CHECK PHASER BASE OSC CONFIG'''
        # check that phaser oscillator amplitude config is valid
        if ((not isinstance(self.ampl_osc_frac_list, list)) or
                (len(self.ampl_osc_frac_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator amplitude array must be list of length {:d}.".format(self._num_phaser_oscs))
        elif sum(self.ampl_osc_frac_list) >= 100.:
            raise ValueError("Error: phaser oscillator amplitudes must sum <100.")

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
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = array([self.freq_phaser_carrier_mhz]) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

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
            raise ValueError("Invalid PSK schedule. All PSK schedules must be of same length.")

        # ensure that sweep targets are lists of appropriate length
        if not (isinstance(self.phase_sweep_arr, list) and (len(self.phase_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid phase_sweep_arr: {:}."
                             "phase_sweep_arr must be list of length {:d}.".format(self.phase_sweep_arr,
                                                                                      self._num_phaser_oscs))
        if not (isinstance(self.freq_sweep_arr, list) and (len(self.freq_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid freq_sweep_arr: {:}."
                             "freq_sweep_arr must be list of length {:d}.".format(self.freq_sweep_arr,
                                                                                      self._num_phaser_oscs))

        # check that waveforms are not too many/not sweeping too hard
        num_waveforms_to_record = len(list(self.phase_osc_sweep_turns_list))
        if num_waveforms_to_record > PULSESHAPER_MAX_WAVEFORMS:
            raise ValueError("Too many waveforms to record ({:d}) - must be fewer than {:d}. "
                             "Reduce length of phase_osc_sweep_turns_list.".format(num_waveforms_to_record,
                                                                                   PULSESHAPER_MAX_WAVEFORMS))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        self.waveform_index_to_compiled_wav = list() # store compiled waveform values
        # note: waveform_index_to_pulseshaper_id is NOT kernel_invariant b/c gets updated in phaser_record
        self.waveform_index_to_pulseshaper_id = zeros(len(self._waveform_param_list), dtype=int32) # store waveform ID linked to DMA sequence

        # calculate block timings and scale amplitudes for ramsey-ing
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        if self.enable_phase_shift_keying:
            # case: PSK w/ramsey delay
            if self.enable_psk_delay:
                num_blocks = 2 * num_psk_blocks - 1
                # note: don't create block_time_list_us here b/c we need to adjust for ramsey delay time sweeps
                block_ampl_scale_list = riffle([1] * num_psk_blocks, [0] * (num_psk_blocks - 1))
            # case: PSK only, no ramsey delay
            else:
                num_blocks = num_psk_blocks
                block_time_list_us = [self.time_heating_us / num_psk_blocks] * num_psk_blocks
                block_ampl_scale_list = [1] * num_psk_blocks
        # case: rabi-style (no PSK, no ramsey delay)
        else:
            num_blocks = 1
            block_time_list_us = [self.time_heating_us]
            block_ampl_scale_list = [1]


        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = array(self.ampl_osc_frac_list)
        _osc_vals_blocks[:, :, 0] *= array([block_ampl_scale_list]).transpose()

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_s_list = array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_sweep_arr * mean(self.freq_osc_sweep_hz_list)) *
                t_update_delay_s_list
        )
        _osc_vals_blocks[:, :, 1] += array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # note: use ::2 since we only update to non-delay blocks
                _osc_vals_blocks[::2, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                      self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                      self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()
            else:
                _osc_vals_blocks[:, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                    self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                    self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()


        '''COMPILE WAVEFORMS SPECIFIC TO EACH PARAMETER'''
        # record phaser waveforms - one for each phase
        for waveform_params in self._waveform_param_list:
            # extract waveform params
            phase_sweep_turns = waveform_params[0]

            # case - PSK + ramsey: create block_time_list_us with target delay time
            if self.enable_phase_shift_keying and self.enable_psk_delay:
                block_time_list_us = riffle([self.time_heating_us / num_psk_blocks] * num_psk_blocks,
                                            [self.time_psk_delay_us] * (num_psk_blocks - 1))

            # create local copy of _osc_vals_blocks and update with target phase
            # note: no need to deep copy b/c it's filled w/immutables
            _osc_vals_blocks_local = np_copy(_osc_vals_blocks)
            _osc_vals_blocks_local[:, :, 1] += self.phase_sweep_arr * phase_sweep_turns

            # specify sequence as a list of blocks, where each block is a dict
            # note: have to instantiate locally each loop b/c dicts aren't deep copied
            _sequence_blocks_local = [
                {
                    "oscillator_parameters": _osc_vals_blocks_local[_idx_block],
                    "config": {
                        "time_us": block_time_list_us[_idx_block],
                        # don't pulse shape for delay blocks lmao
                        "pulse_shaping": self.enable_pulse_shaping and (block_ampl_scale_list[_idx_block] != 0),
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping,
                            "pulse_shape_falling": self.enable_pulse_shaping,
                            "sample_rate_khz": self.freq_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_pulse_shape_rolloff_us
                        }
                    }
                } for _idx_block in range(num_blocks)
            ]
            # create QVSA waveform and store data in a holder
            self.waveform_index_to_compiled_wav.append(
                self.spinecho_wizard.compile_waveform(_sequence_blocks_local)
            )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                8)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # record phaser waveforms
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage during configuration
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.pulse_shaper.waveform_load() # load phaser waveform DMA handles
        _loop_iter = 0  # used to check_termination more frequently

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_sweep_hz =         config_vals[0]
                waveform_num =          int32(config_vals[1])
                phase_ch1_turns =       config_vals[2]
                freq_psrsb_rsb_ftw =    int32(config_vals[3])
                freq_psrsb_carr_ftw =   int32(config_vals[4])
                phas_psrsb_rsb_pow =    int32(config_vals[5])
                phas_psrsb_carr_pow =   int32(config_vals[6])

                # get corresponding waveform parameters and pulseshaper ID from the index
                waveform_params = self._waveform_param_list[waveform_num]
                phaser_waveform = self.waveform_index_to_pulseshaper_id[waveform_num]

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_osc_base_hz_list + freq_sweep_hz * self.freq_sweep_arr
                self.core.break_realtime()
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    self.freq_phaser_carrier_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], 0.],
                    # CH1 phase (via DUC)
                    phase_ch1_turns
                )


                '''STATE PREPARATION'''
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                '''PHASE-SENSITIVE RED SIDEBAND SEQUENCE'''
                # create qvsa displacement
                t_phaser_start_mu = self.phaser_run(phaser_waveform)
                # run PSRSB detection (synchronized to phaser)
                self.psrsb_run(freq_psrsb_rsb_ftw, freq_psrsb_carr_ftw,
                               phas_psrsb_rsb_pow, phas_psrsb_carr_pow, t_phaser_start_mu)

                '''READOUT'''
                self.readout_subsequence.run_dma()
                self.rescue_subsequence.resuscitate()

                '''LOOP CLEANUP'''
                counts = self.readout_subsequence.fetch_count()
                self.rescue_subsequence.detect_death(counts)
                self.update_results(
                    freq_psrsb_rsb_ftw,
                    counts,
                    freq_sweep_hz,
                    waveform_params[0], # note: manually expand waveform_params b/c no variadics in kernel
                    phase_ch1_turns,
                    freq_psrsb_carr_ftw,
                    phas_psrsb_carr_pow,
                    phas_psrsb_carr_pow
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
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TInt64:
        """
        Run the main phaser pulse together with supporting functionality.
        :param waveform_id: the ID of the waveform to run.
        :return:the start time of the phaser oscillator waveform (in machine units, 64b int).
            Useful for phase-coherent device synchronization.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_phaser_mu, self.att_phaser_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        # synchronize to next frame
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
        # record phaser sequences onto DMA for each waveform parameter
        for i in range(len(self._waveform_param_list)):
            # get waveform for given parameters
            # note: use sync RPC to reduce significant overhead of direct data transfer
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self._get_compiled_waveform(i)

            # record phaser pulse sequence and save returned waveform ID
            # note: no need to add slack b/c waveform_record does it for us
            self.waveform_index_to_pulseshaper_id[i] = self.pulse_shaper.waveform_record(
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

    @kernel(flags={"fast-math"})
    def psrsb_run(self, freq_rsb_ftw: TInt32=0, freq_carrier_ftw: TInt32=0,
                  phas_rsb_pow: TInt32=0, phas_carrier_pow: TInt32=0,
                  time_ref_mu: TInt64=-1) -> TNone:
        """
        Run the phase-sensitive red-sideband detection sequence.
        :param freq_rsb_ftw: RSB frequency in FTW.
        :param freq_carrier_ftw: Carrier frequency in FTW.
        :param phas_rsb_pow: RSB phase (relative) in POW.
        :param phas_carrier_pow: Carrier phase (relative) in POW.
        :param time_ref_mu: Fiducial time used to compute coherent/tracking phase updates.
        """
        # set target profile and attenuation
        self.qubit.set_profile(self.profile_729_PSRSB)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_qubit_mu)

        # synchronize start time to coarse RTIO clock
        if time_ref_mu < 0: time_ref_mu = now_mu() & ~0x7

        # run RSB pulse
        self.qubit.set_mu(freq_rsb_ftw, asf=self.ampl_psrsb_rsb_asf,
                          pow_=phas_rsb_pow,
                          profile=self.profile_729_PSRSB,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_ref_mu)
        self.qubit.on()
        delay_mu(self.time_psrsb_rsb_mu)
        self.qubit.off()

        # run carrier pulse
        self.qubit.set_mu(freq_carrier_ftw, asf=self.ampl_psrsb_carr_asf,
                          pow_=phas_carrier_pow,
                          profile=self.profile_729_PSRSB,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_ref_mu)
        self.qubit.on()
        delay_mu(self.time_psrsb_carr_mu)
        self.qubit.off()

