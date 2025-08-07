import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM,
    SidebandReadout, QubitRAP
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes

# todo: add generic all osc use
# todo: change names to make more generic


class SuperDuperResolution(LAXExperiment, Experiment):
    """
    Experiment: Super Duper Resolution

    Supports lots of easily configurable parameter scanning for phaser.
    Experiment name inspired by Sam Crary.
    """
    name = 'Super Duper Resolution'
    kernel_invariants = {
        # hardware values
        'att_eggs_heating_mu', 'freq_superresolution_sweep_hz_list', 'freq_global_offset_hz',
        'freq_superresolution_osc_base_hz_list', 'waveform_index_to_pulseshaper_vals',
        '_waveform_param_list',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP',
        'spinecho_wizard', 'pulse_shaper',

        # RAP
        'att_rap_mu', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',
        '_num_phaser_oscs'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=70, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="Sideband Ratio"))

        self._num_phaser_oscs = 4   # number of phaser oscillators in use
        _argstr = "SDR"  # string to use for argument grouping

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

        # Sideband Readout - extra argument
        self.setattr_argument("time_readout_us_list", Scannable(
                                                        default=[
                                                            ExplicitScan([26.11]),
                                                            RangeScan(0, 800, 200, randomize=True),
                                                        ],
                                                        global_min=1, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group='sideband_readout')

        # RAP-based readout
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=100.7407, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

        # phaser - configurable freq & sweeps
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list", Scannable(
                                                                        default=[
                                                                            ExplicitScan([86.]),
                                                                            CenterScan(86., 0.002, 0.0005, randomize=True),
                                                                        ],
                                                                        global_min=0.005, global_max=4800, global_step=1,
                                                                        unit="MHz", scale=1, precision=6
                                                                    ), group="{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("freq_sweep_arr", PYONValue([-1., 1., 0., 0.]), group="{}.freq_phase_sweep".format(_argstr),
                              tooltip="Defines how oscillator freqs should be adjusted for each value in freq_superresolution_sweep_khz_list."
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the freq value, and osc_1 by -1x the freq value, with the rest untouched."
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("freq_superresolution_sweep_khz_list",    Scannable(
                                                                            default=[
                                                                                ExplicitScan([0.]),
                                                                                CenterScan(0., 4, 0.1, randomize=True),
                                                                                RangeScan(0., 100.0, 26, randomize=True),
                                                                            ],
                                                                            global_min=-10000, global_max=10000, global_step=10,
                                                                            unit="kHz", scale=1, precision=6
                                                                        ), group="{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("phase_sweep_arr", PYONValue([1., 0., 0., 0.]), group="{}.freq_phase_sweep".format(_argstr),
                              tooltip="Defines how oscillator phases should be adjusted for each value in phase_superresolution_sweep_turns_list."
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the phase value, and osc_1 by -1x the phase value, with the rest untouched."
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("phase_superresolution_sweep_turns_list", Scannable(
                                                                            default=[
                                                                                ExplicitScan([0.]),
                                                                                RangeScan(0, 1.0, 26, randomize=True),
                                                                            ],
                                                                            global_min=0.0, global_max=1.0, global_step=1,
                                                                            unit="turns", scale=1, precision=3
                                                                        ), group="{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("phase_eggs_heating_ch1_turns_list",  Scannable(
                                                                        default=[
                                                                            ExplicitScan([0.]),
                                                                            RangeScan(0, 1.0, 26, randomize=True),
                                                                        ],
                                                                        global_min=0.0, global_max=1.0, global_step=1,
                                                                        unit="turns", scale=1, precision=3
                                                                    ), group="{}.freq_phase_sweep".format(_argstr))

        # phaser - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        # phaser - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",  BooleanValue(default=False), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc0_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc1_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc2_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc3_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc4_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("enable_psk_delay", BooleanValue(default=False), group="{}.psk".format(_argstr),
                              tooltip="Add a delay between PSK pulses where oscillator amplitudes are set to 0."
                                      "Can be used to create e.g. a Ramsey or DD-type pulse sequence."
                                      "Requires enable_phase_shift_keying to be enabled; otherwise, does nothing."
                                      "Note: prepare/cleanup methods (e.g. set phaser atts, set ext switch) are not called for the delay.")
        self.setattr_argument("time_psk_delay_us_list", Scannable(
                                                      default=[
                                                          ExplicitScan([26.11]),
                                                          RangeScan(10, 1000, 100, randomize=True),
                                                      ],
                                                      global_min=1, global_max=100000, global_step=1,
                                                      unit="us", scale=1, precision=5
                                                  ), group="{}.psk".format(_argstr),
                              tooltip="Delay time (in us) delay between PSK pulses.")

        # phaser - waveform - general
        self.setattr_argument("time_eggs_heating_us",   NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.waveform".format(_argstr))
        self.setattr_argument("att_eggs_heating_db",    NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.), group="{}.waveform".format(_argstr))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.), group="{}.waveform".format(_argstr))
        self.setattr_argument("freq_superresolution_osc_khz_list",      PYONValue([-702.687, 702.687, 0., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("ampl_superresolution_osc_frac_list",     PYONValue([40., 40., 10., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("phase_superresolution_osc_turns_list",   PYONValue([0., 0., 0., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("phase_oscillators_ch1_offset_turns",     PYONValue([0., 0., 0.5, 0.5, 0.5]), group="{}.waveform".format(_argstr))

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
        self.pulse_shaper = PhaserPulseShaper(self, np.array(self.phase_oscillators_ch1_offset_turns))

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
        # convert values to convenience units
        self.att_eggs_heating_mu = att_to_mu(self.att_eggs_heating_db * dB)
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # format build arguments as numpy arrays with appropriate units
        freq_eggs_carrier_hz_list = np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_superresolution_sweep_hz_list = np.array(list(self.freq_superresolution_sweep_khz_list)) * kHz
        self.freq_superresolution_osc_base_hz_list = np.array(self.freq_superresolution_osc_khz_list) * kHz + self.freq_global_offset_hz
        self.phase_superresolution_sweep_turns_list = np.array(list(self.phase_superresolution_sweep_turns_list))
        self.time_psk_delay_us_list = np.array(list(self.time_psk_delay_us_list))
        self.freq_sweep_arr = np.array(self.freq_sweep_arr, dtype=float)
        self.phase_sweep_arr = np.array(self.phase_sweep_arr, dtype=float)

        # map waveform to index to facilitate waveform recording
        # self.waveform_index_to_phase_sweep_turns = np.arange(len(self.phase_superresolution_sweep_turns_list))
        if self.enable_phase_shift_keying and self.enable_psk_delay:
            waveform_num_list = np.arange(len(self.phase_superresolution_sweep_turns_list) *
                                          len(self.time_psk_delay_us_list))
            self._waveform_param_list = create_experiment_config(self.phase_superresolution_sweep_turns_list,
                                                                 self.time_psk_delay_us_list,
                                                                 shuffle_config=False,
                                                                 config_type=float)
        else:
            waveform_num_list = np.arange(len(self.phase_superresolution_sweep_turns_list))
            self._waveform_param_list = self.phase_superresolution_sweep_turns_list


        '''EXPERIMENT CONFIGURATION'''
        # create config data structure
        # self.config_experiment_list = create_experiment_config(
        #     freq_sideband_readout_ftw_list, freq_eggs_carrier_hz_list, self.freq_superresolution_sweep_hz_list,
        #     self.waveform_index_to_phase_sweep_turns, list(self.phase_eggs_heating_ch1_turns_list),
        #     time_readout_mu_list,
        #     config_type=float, shuffle_config=True
        # )
        self.config_experiment_list = create_experiment_config(
            freq_sideband_readout_ftw_list, freq_eggs_carrier_hz_list, self.freq_superresolution_sweep_hz_list,
            waveform_num_list, list(self.phase_eggs_heating_ch1_turns_list),
            time_readout_mu_list,
            config_type=float, shuffle_config=True
        )

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''CHECK PHASER BASE OSC CONFIG'''
        # check that phaser oscillator amplitude config is valid
        if ((not isinstance(self.ampl_superresolution_osc_frac_list, list)) or
                (len(self.ampl_superresolution_osc_frac_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator amplitude array must be list of length {:d}.".format(self._num_phaser_oscs))
        elif np.sum(self.ampl_superresolution_osc_frac_list) >= 100.:
            raise ValueError("Error: phaser oscillator amplitudes must sum <100.")

        # check that phaser oscillator phase arrays are valid
        if ((not isinstance(self.phase_superresolution_osc_turns_list, list)) or
                (len(self.phase_superresolution_osc_turns_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator phase array must be list of length {:d}.".format(self._num_phaser_oscs))

        # check that phaser oscillator frequencies are valid
        if ((not isinstance(self.freq_superresolution_osc_khz_list, list)) or
                (len(self.freq_superresolution_osc_khz_list) != self._num_phaser_oscs)):
            raise ValueError("Error: phaser oscillator frequency array must be list of length {:d}.".format(self._num_phaser_oscs))
        max_osc_freq_hz = (np.max(list(self.freq_superresolution_sweep_khz_list)) * kHz +
                           max(self.freq_superresolution_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        min_osc_freq_hz = (np.min(list(self.freq_superresolution_sweep_khz_list)) * kHz +
                           max(self.freq_superresolution_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
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
        num_waveforms_to_record = (len(list(self.phase_superresolution_sweep_turns_list)) *
                                   len(list(self.time_psk_delay_us_list)))
        if num_waveforms_to_record > 100:
            raise ValueError("Too many waveforms to record - {:d}. Reduce length of either"
                             "phase_superresolution_sweep_turns_list or time_psk_delay_us_list".format(num_waveforms_to_record))

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        self.waveform_index_to_pulseshaper_vals = list() # store compiled waveform values
        # note: waveform_index_to_pulseshaper_id is NOT kernel_invariant b/c gets updated in phaser_record
        self.waveform_index_to_pulseshaper_id = np.zeros(len(self._waveform_param_list), dtype=np.int32) # store waveform ID linked to DMA sequence

        # calculate block timings and scale amplitudes for ramsey-ing
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        if self.enable_phase_shift_keying:
            # case: PSK w/ramsey delay
            if self.enable_psk_delay:
                num_blocks = 2 * num_psk_blocks - 1
                # block_time_list_us = riffle([self.time_eggs_heating_us / num_psk_blocks] * num_psk_blocks,
                #                             [self.time_psk_delay_us_list[0]] * (num_psk_blocks - 1))
                # note: don't create block_time_list_us here b/c we need to adjust for ramsey delay time sweeps
                block_ampl_scale_list = riffle([1] * num_psk_blocks, [0] * (num_psk_blocks - 1))
            # case: PSK only, no ramsey delay
            else:
                num_blocks = num_psk_blocks
                block_time_list_us = [self.time_eggs_heating_us / num_psk_blocks] * num_psk_blocks
                block_ampl_scale_list = [1] * num_psk_blocks
        # case: rabi-style (no PSK, no ramsey delay)
        else:
            num_blocks = 1
            block_time_list_us = [self.time_eggs_heating_us]
            block_ampl_scale_list = [1]


        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = np.zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = np.array(self.ampl_superresolution_osc_frac_list)
        _osc_vals_blocks[:, :, 0] *= np.array([block_ampl_scale_list]).transpose()

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_s_list = np.array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        phase_osc_update_delay_turns_list = (
                (self.freq_superresolution_osc_base_hz_list +
                 self.freq_sweep_arr * np.mean(self.freq_superresolution_sweep_hz_list)) *
                t_update_delay_s_list # account for relative delay period between each osc
        )
        _osc_vals_blocks[:, :, 1] += (np.array(self.phase_superresolution_osc_turns_list) +
                                      phase_osc_update_delay_turns_list)

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # note: use ::2 since we only update to non-delay blocks
                _osc_vals_blocks[::2, :, 1] += np.array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                         self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                         self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()
            else:
                _osc_vals_blocks[:, :, 1] += np.array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                       self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                       self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()


        '''COMPILE WAVEFORMS SPECIFIC TO EACH PARAMETER'''
        # record phaser waveforms - one for each phase
        # for phase in self.phase_superresolution_sweep_turns_list:
        for waveform_params in self._waveform_param_list:
            # extract waveform params
            phase_sweep_turns = waveform_params[0]
            time_us_delay = waveform_params[1]

            # case - PSK + ramsey: create block_time_list_us with target delay time
            if self.enable_phase_shift_keying and self.enable_psk_delay:
                block_time_list_us = riffle([self.time_eggs_heating_us / num_psk_blocks] * num_psk_blocks,
                                            [time_us_delay] * (num_psk_blocks - 1))

            # create local copy of _osc_vals_blocks and update with target phase
            # note: no need to deep copy b/c it's filled w/immutables
            _osc_vals_blocks_local = np.copy(_osc_vals_blocks)
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
            self.waveform_index_to_pulseshaper_vals.append(
                self.spinecho_wizard.compile_waveform(_sequence_blocks_local)
            )

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
                # phase_sweep_idx =   np.int32(config_vals[3])
                waveform_idx =      np.int32(config_vals[3])
                phase_ch1_turns =   config_vals[4]
                time_readout_mu =   np.int64(config_vals[5])

                # get corresponding waveform and pulseshaper ID from the index
                # phase_sweep_turns = self.phase_superresolution_sweep_turns_list[phase_sweep_idx]
                waveform_params = self._waveform_param_list[waveform_idx]
                phaser_waveform = self.waveform_index_to_pulseshaper_id[waveform_idx]

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_superresolution_osc_base_hz_list + freq_sweep_hz * self.freq_sweep_arr
                self.core.break_realtime()
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    carrier_freq_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], 0.],
                    phase_ch1_turns
                )
                delay_mu(35000)

                # set qubit readout frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_sb_readout, phase_mode=PHASE_MODE_CONTINUOUS)

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                '''QVSA PULSE'''
                self.phaser_run(phaser_waveform)

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
                    # phase_sweep_turns,
                    waveform_params[0],
                    waveform_params[1],
                    phase_ch1_turns,
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
        Run the main EGGS pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_eggs_heating_mu, self.att_eggs_heating_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # EGGS - STOP
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each waveform parameter
        # for i in range(len(self.phase_superresolution_sweep_turns_list)):
        for i in range(len(self._waveform_param_list)):
            # get waveform for given sweep phase
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self.waveform_index_to_pulseshaper_vals[i]

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            self.waveform_index_to_pulseshaper_id[i] = self.pulse_shaper.waveform_record(
                _wav_data_ampl, _wav_data_phas, _wav_data_time
            )
            self.core.break_realtime()

