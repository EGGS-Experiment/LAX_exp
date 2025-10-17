from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_TRACKING

from numpy import int32, int64, array, zeros

from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, Readout, SidebandCoolContinuousRAM, FockOverlap

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

# todo: add cat ampl


class ContinuousSamplingAmplRDX(LAXExperiment, Experiment):
    """
    Experiment: Continuous Sampling Ampl RDX

    Synchronized/correlation spectroscopy using a dynamical-decoupling protocol w/QVSA
        for arbitrary frequency sensing.
    Uses burst sequencing/readout to reduce latency/improve sample rate.
    Uses fock state-based quantum amplification (following https://www.nature.com/articles/s41467-019-10576-4)
        with only RAP-based primitives (to mitigate pulse errors).
    """
    name = 'Continuous Sampling Ampl'
    kernel_invariants = {
        # hardware values
        'freq_osc_base_hz_list', 'freq_phaser_carrier_hz', 'att_phaser_mu', 'pulseshaper_vals',
        'sample_period_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'fock_subsequence', 'spinecho_wizard', 'pulse_shaper',

        # configs
        'profile_729_SBC', 'profile_729_RAP', 'profile_729_fock',
        '_num_phaser_oscs', '_enable_osc_clr', '_burst_samples', '_num_bursts', '_idx_burst_samples',
        'TIME_SAMPLE_PERIOD_MIN_SLACK_MU',

        # PSRSB extra
        'protocol_type_num', 'profile_729_PSRSB', 'att_psrsb_mu',
        'ampl_psrsb_rsb_asf', 'ampl_psrsb_carr_asf', 'time_psrsb_rsb_mu', 'time_psrsb_carr_mu',
        'freq_psrsb_rsb_ftw', 'freq_psrsb_carr_ftw', 'phas_psrsb_rsb_pow', 'phas_psrsb_carr_pow',
    }

    def build_experiment(self):
        # exp-specific variables
        self._num_phaser_oscs = 5   # num phaser oscillators in use
        self._burst_samples =   50  # num experiment shots in a "burst" - must be <100 b/c EdgeCounter limits
        self.TIME_SAMPLE_PERIOD_MIN_SLACK_MU = 125000   # min slack between actual shot time and sample period to prevent RTIOUnderflows

        # core arguments
        self.setattr_argument("num_samples",        NumberValue(default=2000, precision=0, step=1, min=1, max=10000000),
                              tooltip="Number of samples to record.")
        self.setattr_argument("sample_period_ms",   NumberValue(default=11., precision=6, min=5, max=1e5, step=1, unit="ms", scale=1.),
                              tooltip="Overall time for each sample."
                                      "Must be greater than the actual pulse time by some minimum slack.")
        self.setattr_argument("protocol_type", EnumerationValue(["FOCK", "PSRSB"], default="FOCK"))

        # allocate relevant beam profiles
        self.profile_729_SBC =      1
        self.profile_729_RAP =      2
        self.profile_729_fock =     3
        self.profile_729_PSRSB =    4

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        # quantum amplification: fock overlap
        # note: FockOverlap handles the RAP pulse config for us
        self.fock_subsequence = FockOverlap(
            self, ram_profile=self.profile_729_RAP,
            ram_addr_start=202, num_samples=250, pulse_shape="blackman"
        )

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

        # set build arguments
        self._build_arguments_PSRSB()
        self._build_arguments_waveform()
        self._build_arguments_modulation()

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def _build_arguments_PSRSB(self):
        """
        Set specific arguments for PSRSB beam configuration.
        """
        _argstr = "PSRSB"    # string to use for arguments

        # PSRSB - general
        self.setattr_argument("att_psrsb_db", NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group='{}'.format(_argstr))

        # PSRSB - RSB pulse
        self.setattr_argument("ampl_psrsb_rsb_pct", NumberValue(default=50, precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group='{}'.format(_argstr))
        self.setattr_argument("time_psrsb_rsb_us",  NumberValue(default=25.41, precision=3, step=5, min=1, max=10000000, scale=1., unit="us"),
                              group='{}'.format(_argstr))
        self.setattr_argument("freq_psrsb_rsb_mhz", NumberValue(default=100.7594, precision=6, step=5, min=60., max=200., scale=1., unit="MHz"),
                              group='{}'.format(_argstr))
        self.setattr_argument("phas_psrsb_rsb_turns",   NumberValue(default=0., precision=5, step=0.05, min=-1.1, max=1.1, scale=1., unit="turns"),
                              group='{}'.format(_argstr))

        # PSRSB - carrier pulse
        self.setattr_argument("ampl_psrsb_carr_pct",    NumberValue(default=50, precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group='{}'.format(_argstr))
        self.setattr_argument("time_psrsb_carr_us",     NumberValue(default=1.21, precision=3, step=1, min=1, max=10000000, scale=1., unit="us"),
                              group='{}'.format(_argstr))
        self.setattr_argument("freq_psrsb_carr_mhz",    NumberValue(default=101.0881, precision=6, step=5, min=60., max=200., scale=1., unit="MHz"),
                              group='{}'.format(_argstr))
        self.setattr_argument("phas_psrsb_carr_turns",  NumberValue(default=0., precision=5, step=0.05, min=-1.1, max=1.1, scale=1., unit="turns"),
                              group='{}'.format(_argstr))

    def _build_arguments_waveform(self):
        """
        Build core waveform arguments for the QVSA pulse.
        """
        _argstr = "CS"  # string to use for arguments

        # waveform - global config
        self.setattr_argument("att_phaser_db",
                              NumberValue(default=28., precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Phaser attenuation to be used for both CH0 and CH1.")
        self.setattr_argument("freq_phaser_carrier_mhz",
                              NumberValue(default=86., precision=7, step=1, min=0.001, max=4800, unit="MHz", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Phaser output center frequency.\n"
                                      "Note: actual center frequency depends on the devices.phaser.freq_center_mhz dataset argument, "
                                      "which should be manually entered into the dataset manager by the user after "
                                      "configuring the TRF and NCO via e.g. the phaser_configure tool.\n"
                                      "Ensure all values are set correctly.")
        self.setattr_argument("freq_global_offset_mhz",
                              NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Apply a frequency offset via the phaser oscillators to avoid any DUC/NCO/TRF output spurs.\n"
                                      "Range is limited by the phaser oscillator freq range, i.e. [-10, 10] MHz (includes the frequencies in freq_osc_khz_list).")
        self.setattr_argument("phase_global_ch1_turns",
                              NumberValue(default=0., precision=5, step=0.05, min=-1.1, max=1.1, unit="turns", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Sets a global CH1 phase via the DUC.\n"
                                      "Note: the eggs.phas_ch1_inherent_turns dataset argument is overridden "
                                      "in this experiment.")
        self.setattr_argument("osc_num_target_list", PYONValue([]),
                              group="{}.global".format(_argstr),
                              tooltip="Select oscillators which should reset their phase at the beginning of each shot.\n"
                                      "All oscillators excluded from this list are continuously phase-tracking by default.")

        # waveform - custom specification
        self.setattr_argument("time_osc_pulse_us",
                              NumberValue(default=500, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Time for a SINGLE SEGMENT of the pulse."
                                      "e.g. a time_osc_pulse_us of 500us with 2 segments => total time of 1ms.")
        self.setattr_argument("freq_osc_khz_list", PYONValue([-702.687, 702.687, 0.005, 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator frequencies.")
        self.setattr_argument("ampl_osc_frac_list", PYONValue([25., 25., 49., 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator amplitudes. Applied to both CH0 and CH1.\n"
                                      "Note: CH1 amplitudes will be scaled by the amplitude scaling factors in devices.phaser.ch1.ampl_ch1_osc_scale_arr.")
        self.setattr_argument("phase_osc_turns_list", PYONValue([0., 0., 0., 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Relative phases between each phaser oscillator. Applied on both CH0 and CH1.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0.5, 0.5, 0.5]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Sets the relative CH1 phase via the phaser oscillators.")

    def _build_arguments_modulation(self):
        """
        Build core modulation arguments for the QVSA pulse.
        """
        _argstr = "CS"  # string to use for arguments

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=True),
                              group='{}.shape'.format(_argstr),
                              tooltip="Applies pulse shaping to the edges of the phaser pulse.\n"
                                      "Note: pulse shaping is applied to each constituent PSK block.")
        self.setattr_argument("type_pulse_shape",   EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.shape'.format(_argstr),
                              tooltip="Pulse shape type to be used.")
        self.setattr_argument("time_pulse_shape_rolloff_us",
                              NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Time constant of the pulse shape. This is used for both the pulse rollon AND rolloff.\n"
                                      "e.g. a 1ms main pulse time with 100us time_pulse_shape_rolloff_us will result in a 1ms + 2*100us = 1.2ms total pulse time.\n"
                                      "All constituent PSK blocks will have this pulse time applied.\n"
                                      "Note: DMA issues limit the total number of samples (i.e. time_pulse_shape_rolloff_us * freq_pulse_shape_sample_khz).")
        self.setattr_argument("freq_pulse_shape_sample_khz",
                              NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Sample rate used for pulse shaping.\n"
                                      "This value is inexact and is fixed at multiples of the phaser oscillator update "
                                      "rate (i.e. 40ns) times the number of oscillators in use.")

        # waveform - PSK (Phase-shift Keying) - general
        self.setattr_argument("enable_phase_shift_keying", BooleanValue(default=False),
                              group="{}.psk".format(_argstr),
                              tooltip="Enable PSK-ing: break the main pulse into individual blocks with different phases.\n"
                                      "Number of PSKs is determined by number of phases in phase_osc<x>_psk_turns. "
                                      "All oscillator PSK schedules must have same length.")
        self.setattr_argument("enable_psk_delay", BooleanValue(default=False),
                              group="{}.psk".format(_argstr),
                              tooltip="Add a delay between PSK pulses where oscillator amplitudes are set to 0. "
                                      "Can be used to create e.g. a Ramsey or DD-type pulse sequence.\n"
                                      "Requires enable_phase_shift_keying to be enabled; otherwise, does nothing.\n"
                                      "Note: prepare/cleanup methods (e.g. set phaser atts, set ext switch) are not called for the delay.")
        self.setattr_argument("time_psk_delay_us",
                              NumberValue(default=200., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.),
                              group="{}.psk".format(_argstr),
                              tooltip="Delay time (in us) delay between PSK pulses. Used for ramsey-ing.\n"
                                      "Note: enable_phase_shift_keying AND enable_psk_delay must be True.")

        # waveform - PSK (Phase-shift Keying) - schedule
        self.setattr_argument("phase_osc0_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc0.")
        self.setattr_argument("phase_osc1_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc1.")
        self.setattr_argument("phase_osc2_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc2.")
        self.setattr_argument("phase_osc3_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc3.")
        self.setattr_argument("phase_osc4_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr), tooltip="PSK phase schedule for osc4.")

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''
        GENERAL SETUP
        '''
        self._prepare_argument_checks() # check input arguments for safety

        # process protocol type (b/c don't want to compare str in kernel)
        if self.protocol_type == "FOCK":    self.protocol_type_num = 0
        elif self.protocol_type == "PSRSB": self.protocol_type_num = 1

        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks stored as [block_num, osc_num] and store [ampl_pct, phase_turns]
        #   e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_osc_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, array(self.phase_osc_ch1_offset_turns))


        '''
        BURST DMA PREPARATION
        '''
        ### for convenience/speed, rigidly group shots into a "burst" ###
        self._num_bursts = round(self.num_samples / self._burst_samples)    # for convenience, ensure num_samples is integer multiple of bursts
        self.num_samples = self._num_bursts * self._burst_samples # NOTE: REDEFINITION of the num_samples argument
        
        # protected variables for coredevice use
        self._counts_burst = zeros(self._burst_samples, dtype=int32)      # store burst counts
        self._times_start_burst = zeros(self._burst_samples, dtype=int64) # store burst start timestamps
        self._times_stop_burst = zeros(self._burst_samples, dtype=int64)  # store burst stop timestamps
        self._sequence_dma_handle = (0, int64(0), int32(0), False)        # store sequence DMA handle
        
        self._t_start_mu = int64(0) # store timestamp of start of each shot
        self._idx_burst_samples = list(range(self._burst_samples)) # hold iterators for each burst loop (hopefully reduces overhead)
        self._time_exp_shot_mu = int64(0)   # measure actual shot period (during DMA recording) to ensure user-specified sample period is viable


        '''
        HARDWARE VALUES - READOUT
        '''
        # PSRSB DDS waveform params
        self.att_psrsb_mu = self.qubit.cpld.att_to_mu(self.att_psrsb_db * dB)

        self.ampl_psrsb_rsb_asf =   self.qubit.amplitude_to_asf(self.ampl_psrsb_rsb_pct / 100.)
        self.ampl_psrsb_carr_asf =  self.qubit.amplitude_to_asf(self.ampl_psrsb_carr_pct / 100.)

        self.time_psrsb_rsb_mu =    self.core.seconds_to_mu(self.time_psrsb_rsb_us * us)
        self.time_psrsb_carr_mu =   self.core.seconds_to_mu(self.time_psrsb_carr_us * us)

        self.freq_psrsb_rsb_ftw =   self.qubit.frequency_to_ftw(self.freq_psrsb_rsb_mhz * MHz)
        self.freq_psrsb_carr_ftw =  self.qubit.frequency_to_ftw(self.freq_psrsb_carr_mhz * MHz)

        self.phas_psrsb_rsb_pow =   self.qubit.turns_to_pow(self.phas_psrsb_rsb_turns)
        self.phas_psrsb_carr_pow =  self.qubit.turns_to_pow(self.phas_psrsb_carr_turns)


        '''
        HARDWARE VALUES - CONFIG
        '''
        # ensure sample interval is multiple of phaser frame period (320ns)
        self.sample_period_mu = int64(
            round(
                self.core.seconds_to_mu(self.sample_period_ms * ms) /
                self.phaser_eggs.t_frame_mu
            ) * self.phaser_eggs.t_frame_mu
        )

        # convert build arguments to appropriate formats
        self.att_phaser_mu = att_to_mu(self.att_phaser_db * dB)
        self.freq_phaser_carrier_hz = (self.freq_phaser_carrier_mhz - self.freq_global_offset_mhz) * MHz
        self.freq_osc_base_hz_list = array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_mhz * MHz

        self._prepare_waveform()    # configure waveform via pulse shaper & spin echo wizard

        # account for case when osc_num_target_list is empty
        if len(self.osc_num_target_list) == 0:
            self.osc_num_target_list = [-1]
            self._enable_osc_clr = False
        else:   self._enable_osc_clr = True

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''
        CHECK PHASER BASE OSC CONFIG
        '''
        # check that phaser oscillator amplitude config is valid
        if (not isinstance(self.ampl_osc_frac_list, list)) or (len(self.ampl_osc_frac_list) != self._num_phaser_oscs):
            raise ValueError("Invalid argument: phaser osc amplitude array must be list of length {:d}.".format(
                self._num_phaser_oscs))
        elif sum(self.ampl_osc_frac_list) >= 100.:
            raise ValueError("Invalid argument: phaser osc amplitudes must sum <100.")

        # check that phaser oscillator phase arrays are valid
        if (not isinstance(self.phase_osc_turns_list, list)) or (len(self.phase_osc_turns_list) != self._num_phaser_oscs):
            raise ValueError("Invalid argument: phaser osc phase array must be list of length {:d}.".format(
                self._num_phaser_oscs))

        # check that phaser oscillator frequencies are valid
        if (not isinstance(self.freq_osc_khz_list, list)) or (len(self.freq_osc_khz_list) != self._num_phaser_oscs):
            raise ValueError("Invalid argument: phaser osc frequency array must be list of length {:d}.".format(
                self._num_phaser_oscs))
        max_osc_freq_hz = max(self.freq_osc_khz_list) * kHz + (self.freq_global_offset_mhz * MHz)
        min_osc_freq_hz = min(self.freq_osc_khz_list) * kHz + (self.freq_global_offset_mhz * MHz)
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Invalid argument: phaser osc frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_carrier_freq_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_phaser_carrier_mhz * MHz)
        if phaser_carrier_freq_dev_hz >= 300. * MHz:
            raise ValueError("Invalid argument: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")


        '''
        CHECK PHASER WAVEFORM CONFIG
        '''
        # check that PSK schedule is valid
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (not isinstance(psk_schedule, list)) or (len(psk_schedule) != num_psk_blocks)
            for psk_schedule in (
                self.phase_osc0_psk_turns, self.phase_osc1_psk_turns, self.phase_osc2_psk_turns,
                self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule. All PSK schedules must be of same length.")

        # ensure that spinecho-ing makes sense
        if self.enable_psk_delay and not self.enable_phase_shift_keying:
            raise ValueError("Invalid waveform configuration. Cannot have delays enabled without PSKing.")

        # ensure osc_num_target_list is a valid list
        if not (isinstance(self.osc_num_target_list, list)):
            raise ValueError("Invalid osc_num_target_list. Must be of type list.")
        # ensure that osc_num_target_list contains a valid selection of oscillators
        elif (len(self.osc_num_target_list) != 0) and not all((
                all(isinstance(val, int) for val in self.osc_num_target_list),
                max(self.osc_num_target_list) <= self._num_phaser_oscs,
                min(self.osc_num_target_list) >= 0,
                len(self.osc_num_target_list) <= self._num_phaser_oscs,
                len(set(self.osc_num_target_list)) <= self._num_phaser_oscs
        )):
            raise ValueError("Invalid osc_num_target_list. "
                             "Must be a list of fewer than {:d} numbers in [0, {:d}].".format(
                self._num_phaser_oscs, self._num_phaser_oscs-1))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizardRDX and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.pulseshaper_vals = None  # store compiled waveforms from pulseshaper
        self.pulseshaper_id = int32(0)  # store waveform ID for pulseshaper

        # calculate block timings and scale amplitudes for ramsey-ing
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                num_blocks = 2 * num_psk_blocks - 1
                block_time_list_us = riffle([self.time_osc_pulse_us] * num_psk_blocks,
                                            [self.time_psk_delay_us] * (num_psk_blocks - 1))
                block_ampl_scale_list = riffle([1] * num_psk_blocks, [0] * (num_psk_blocks - 1))
            else:
                num_blocks = num_psk_blocks
                block_time_list_us = [self.time_osc_pulse_us] * num_psk_blocks
                block_ampl_scale_list = [1] * num_psk_blocks
        else:
            num_blocks = 1
            block_time_list_us = [self.time_osc_pulse_us]
            block_ampl_scale_list = [1]

        '''PROGRAM & COMPILE WAVEFORM'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = array(self.ampl_osc_frac_list)
        _osc_vals_blocks[:, :, 0] *= array([block_ampl_scale_list]).transpose()

        # set oscillator phases and account for oscillator update delays
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_s_list = array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        _osc_vals_blocks[:, :, 1] += (array(self.phase_osc_turns_list) +
                                      self.freq_osc_base_hz_list * t_update_delay_s_list)

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # note: use ::2 since we only update to non-delay blocks
                _osc_vals_blocks[::2, :, 1] += array([
                                                            self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                            self.phase_osc2_psk_turns,
                                                            self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
                                                        ][:self._num_phaser_oscs]).transpose()
            else:
                _osc_vals_blocks[:, :, 1] += array([
                                                          self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                          self.phase_osc2_psk_turns,
                                                          self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
                                                      ][:self._num_phaser_oscs]).transpose()

        # specify sequence as a list of blocks, where each block is a dict
        _sequence_blocks = [
            {
                "oscillator_parameters": _osc_vals_blocks[_idx_block],
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
        self.pulseshaper_vals = self.spinecho_wizard.compile_waveform(_sequence_blocks)

    @property
    def results_shape(self):
        return (self.num_samples, 3)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        """
        Initialize relevant components of the experiment in-kernel immediately before run.
        """
        # get phaser waveform for PulseShaper
        # note: we don't do any error checking here, so have to be really sure all OK
        ampl_frac_list, phas_turns_list, sample_interval_mu_list = self.pulseshaper_vals
        self.core.break_realtime()

        # record full sequence for burst reasons
        delay_mu(1000000)  # add slack - 1ms
        # align to phaser frame for later deterministic playback
        # note: only need to align to frame b/c frame is multiple of coarse RTIO
        at_mu(self.phaser_eggs.get_next_frame_mu())
        with self.core_dma.record('_SEQUENCE_SHOT'):

            '''INITIALIZE ION'''
            t0 = now_mu()   # record exp shot period - start time
            self.initialize_subsequence.run()   # doppler cool & prep spin into |S-1/2, mJ=-1/2>
            self.sidebandcool_subsequence.run() # multimode cont. SBC w/854nm quenching

            '''PREPARE MOTIONAL STATE/AMPLIFICATION'''
            if self.protocol_type_num == 0:
                self.fock_subsequence.run_fock_generate() # generate higher fock state

            '''WAVEFORM SEQUENCE'''
            # prepare phaser for output
            self.phaser_eggs.phaser_setup(self.att_phaser_mu, self.att_phaser_mu)
            t_phaser_start_mu = self.phaser_eggs.get_next_frame_mu()   # get fiducial time for PSRSB

            # run oscillator waveform
            at_mu(t_phaser_start_mu)
            for i in range(len(ampl_frac_list)):
                self.pulse_shaper._waveform_point(ampl_frac_list[i], phas_turns_list[i])
                delay_mu(sample_interval_mu_list[i])
            # stop phaser output
            self._phaser_stop_rdx()

            '''READOUT'''
            # higher fock overlap readout
            if self.protocol_type_num == 0:
                self.fock_subsequence.run_fock_read()
            # PSRSB readout
            elif self.protocol_type_num == 2:
                self.psrsb_run(t_phaser_start_mu)

            self.readout_subsequence.run() # state-selective fluorescence readout

            '''CLEAN UP - SET RESCUE UNTIL NEXT SHOT'''
            # tmp remove
            self.pump.rescue()
            self.pump.on()
            self.repump_cooling.on()
            self.repump_qubit.on()
            # tmp remove

            t1 = now_mu()   # record exp shot period - stop time

        # record exp shot period
        self._time_exp_shot_mu = t1 - t0
        print("\n\t\tContSamp shot period (ms):  ", self.core.mu_to_seconds(t1 - t0) / ms, "\n")
        self.core.break_realtime()


        ### PHASER INITIALIZATION ###
        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_att_mu(0x00)

        # configure global phaser configs (e.g. DUC)
        self.phaser_eggs.frequency_configure(
            self.freq_phaser_carrier_hz,  # carrier frequency (via DUC)
            # oscillator frequencies
            [self.freq_osc_base_hz_list[0], self.freq_osc_base_hz_list[1], self.freq_osc_base_hz_list[2],
             self.freq_osc_base_hz_list[3], self.freq_osc_base_hz_list[4]],
            self.phase_global_ch1_turns  # global CH1 phase
        )

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # ensure shot period is valid and that we have sufficient slack between shots
        # note: do in run_main instead of initialize_experiment to ensure that cleanup happens
        if self.sample_period_mu - self._time_exp_shot_mu < self.TIME_SAMPLE_PERIOD_MIN_SLACK_MU: # 125us (i.e. a break_realtime)
            print("\n\t\tError: actual shot period exceeds allotted sample period.\n")
            self.core.break_realtime()
            return
        delay_mu(125000)

        # load burst sequence DMA handle
        self._sequence_dma_handle = self.core_dma.get_handle('_SEQUENCE_SHOT')
        self.core.break_realtime()  # add slack

        # other setup
        delay_mu(10000000) # add 10ms slack before beginning (just in case)
        self._t_start_mu = self.phaser_eggs.get_next_frame_mu()  # ensure samples are evenly spaced & aligned to phaser frame

        # MAIN LOOP
        for burst_num in range(self._num_bursts):
            # burst sequence for low latency
            self.burst_sequence()
            # burst readout & store results quickly
            self.burst_readout()

            # # check termination more frequently in case reps are low
            # self.check_termination()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def burst_sequence(self) -> TNone:
        """
        Submit a number of experimental sequences in "burst" format.
        """
        # run a number of shots in a "burst" (for latency/slack)
        for sample_num in self._idx_burst_samples:
            self._t_start_mu += self.sample_period_mu
            # note: no need to align to coarse RTIO/phase frame b/c sample_period_mu
            #   guaranteed to be multiple of phaser frame period
            at_mu(self._t_start_mu)
            self.core_dma.playback_handle(self._sequence_dma_handle)

            # record start and stop times
            self._times_start_burst[sample_num] = self._t_start_mu
            self._times_stop_burst[sample_num] = now_mu()

    @kernel(flags={"fast-math"})
    def burst_readout(self) -> TNone:
        """
        Read a number of PMT counts in "burst" format.
        Values are stored in self._burst_samples.
        """
        # burst readout
        for sample_num in self._idx_burst_samples:
            self._counts_burst[sample_num] = self.readout_subsequence.fetch_count()
        # note: add slack all at once instead of inside loop to reduce overhead
        delay_mu(50000)

        # update dataset
        self.update_results(self._times_start_burst,
                            self._counts_burst,
                            self._times_stop_burst)

    @rpc(flags={"async"})
    def update_results(self, time_start_arr_mu: TArray(TInt64), count_arr: TArray(TInt32),
                       time_stop_arr_mu: TArray(TInt64)) -> TNone:
        """
        Overload LAX_exp.base_experiment to allow burst updates.
        Records data from the main sequence in the experiment dataset.
        :param time_start_arr_mu: array of start times (in machine units) for each shot
        :param count_arr: array of measured counts for each shot
        :param time_stop_arr_mu: array of start times (in machine units) for each shot
        """
        # store results in main dataset
        res_arr = array([time_start_arr_mu, count_arr, time_stop_arr_mu]).transpose()
        num_values = len(res_arr)
        self.mutate_dataset('results', (self._result_iter, self._result_iter + num_values), res_arr)

        # get select values for count monitoring
        counts_tmp = array([
            count_val
            for idx, count_val in enumerate(count_arr)
            if (idx + self._result_iter) % self._dynamic_reduction_factor == 0
        ])
        num_counts = len(counts_tmp)
        self.mutate_dataset('temp.counts.trace',
                            (self._counts_iter, self._counts_iter + num_counts),
                            counts_tmp)
        self._counts_iter += num_counts

        # update completion status
        self.set_dataset('management.dynamic.completion_pct',
                         round(self._result_iter * self._completion_iter_to_pct, 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += num_values

    @kernel(flags={"fast-math"})
    def _phaser_stop_rdx(self) -> TNone:
        """
        Stop all output on phaser WITHOUT clearing certain phase accumulators.
        Custom function is needed b/c phaser_eggs.phaser_stop() DOES clear the phase accumulators by default.
        Sets maximum attenuation to prevent output leakage and pulses relevant TTLs with the appropriate delays.
        """
        # set amplitudes of all oscillators to 0
        for i in range(5):
            with parallel:
                self.phaser_eggs.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=0)
                self.phaser_eggs.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=0)
                delay_mu(self.phaser_eggs.t_sample_mu)

        # clear only phase accumulators of target oscs
        if self._enable_osc_clr:
            for osc_num in self.osc_num_target_list:
                with parallel:
                    self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(amplitude=0., phase=0., clr=1)
                    self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(amplitude=0., phase=0., clr=1)
                    delay_mu(self.phaser_eggs.t_sample_mu)

        # add delay for oscillator updates to account for pipeline latency
        delay_mu(2560)  # 8 frame periods - 2.56us
        # stop phaser amp switches & deactivate integrator hold
        with parallel:
            self.phaser_eggs.ch0_amp_sw.off()
            self.phaser_eggs.ch1_amp_sw.off()
            self.phaser_eggs.int_hold.off()

        # add delay time after EGGS pulse to allow RF servo to re-lock
        delay_mu(self.phaser_eggs.time_phaser_holdoff_mu)

        # switch off EGGS attenuators to prevent phaser leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def psrsb_run(self, time_ref_mu: TInt64=-1) -> TNone:
        """
        Run the phase-sensitive red-sideband detection sequence.
        :param time_ref_mu: Fiducial time used to compute coherent/tracking phase updates.
        """
        # set target profile and attenuation
        self.qubit.set_profile(self.profile_729_PSRSB)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_psrsb_mu)

        # synchronize start time to coarse RTIO clock
        if time_ref_mu < 0: time_ref_mu = now_mu() & ~0x7

        # run RSB pulse
        self.qubit.set_mu(self.freq_psrsb_rsb_ftw, asf=self.ampl_psrsb_rsb_asf,
                          pow_=self.phas_psrsb_rsb_pow,
                          profile=self.profile_729_PSRSB,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_ref_mu)
        self.qubit.on()
        delay_mu(self.time_psrsb_rsb_mu)
        self.qubit.off()

        # run carrier pulse
        self.qubit.set_mu(self.freq_psrsb_carr_ftw, asf=self.ampl_psrsb_carr_asf,
                          pow_=self.phas_psrsb_carr_pow,
                          profile=self.profile_729_PSRSB,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_ref_mu)
        self.qubit.on()
        delay_mu(self.time_psrsb_carr_mu)
        self.qubit.off()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        Analyze experimental results.
        """
        # save actual shot period
        self.set_dataset("exp_shot_period_ms", self.core.mu_to_seconds(self._time_exp_shot_mu) / ms)
        # todo: fft stuff etc.
        # todo: report stats - peak/bgr, variability
        pass

