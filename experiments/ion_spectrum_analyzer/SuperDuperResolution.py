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

# todo: add configurable delay for waveform - allow absorption of ramsey-ing
# todo: add generic all osc use
# todo: make freq sweep a specifiable array instead of fixed values
# todo: update argument checks to ensure comprehensive - ensure types all OK
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
        'att_eggs_heating_mu', 'freq_superresolution_sweep_hz_list', 'freq_update_arr', 'freq_global_offset_hz',
        'phase_superresolution_sweep_turns_list', 'freq_superresolution_osc_base_hz_list',

        # EGGS/phaser related
        'waveform_index_to_pulseshaper_vals', 'waveform_index_to_phase_sweep_turns',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP', 'spinecho_wizard',

        # configs
        'profile_729_sb_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list',
        '_num_phaser_oscs'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=70, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type",   EnumerationValue(["Sideband Ratio", "RAP"], default="Sideband Ratio"))
        self.setattr_argument("randomize_config", BooleanValue(default=True))

        self._num_phaser_oscs = 4   # number of phaser oscillators in use
        _argstr = "SDR"  # string to use for arguments

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
        self.setattr_argument("att_rap_db", NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct", NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz", NumberValue(default=100.7407, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz", NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us", NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

        # superresolution - configurable freq & sweeps
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list", Scannable(
                                                                        default=[
                                                                            ExplicitScan([86.]),
                                                                            CenterScan(86., 0.002, 0.0005, randomize=True),
                                                                        ],
                                                                        global_min=0.005, global_max=4800, global_step=1,
                                                                        unit="MHz", scale=1, precision=6
                                                                    ), group="{}.freq_phase_sweep".format(_argstr))

        self.setattr_argument("target_freq_sweep",  EnumerationValue(['Secular', 'Probe', 'X1', 'X2'],
                                                                     default='Secular'), group = "{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("freq_superresolution_sweep_khz_list",    Scannable(
                                                                            default=[
                                                                                ExplicitScan([0.]),
                                                                                CenterScan(0., 4, 0.5, randomize=True),
                                                                            ],
                                                                            global_min=-10000, global_max=10000, global_step=10,
                                                                            unit="kHz", scale=1, precision=6
                                                                        ), group = "{}.freq_phase_sweep".format(_argstr))

        self.setattr_argument("target_phase_sweep", EnumerationValue(['RSB', 'BSB', 'RSB - BSB',
                                                                      'Carrier 0', 'Carrier 1', 'Carrier 0 - Carrier 1'],
                                                                     default='RSB'), group = "{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("phase_superresolution_sweep_turns_list", Scannable(
                                                                            default=[
                                                                                ExplicitScan([0.]),
                                                                                RangeScan(0, 1.0, 3, randomize=True),
                                                                            ],
                                                                            global_min=0.0, global_max=1.0, global_step=1,
                                                                            unit="turns", scale=1, precision=3
                                                                        ), group = "{}.freq_phase_sweep".format(_argstr))
        self.setattr_argument("phase_eggs_heating_ch1_turns_list",  Scannable(
                                                                        default=[
                                                                            ExplicitScan([0.]),
                                                                            RangeScan(0, 1.0, 21, randomize=True),
                                                                        ],
                                                                        global_min=0.0, global_max=1.0, global_step=1,
                                                                        unit="turns", scale=1, precision=3
                                                                    ), group = "{}.freq_phase_sweep".format(_argstr))

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1000, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        # EGGS RF - waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",  BooleanValue(default=False), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_superresolution_rsb_psk_turns",    PYONValue([0., 0.5]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_superresolution_bsb_psk_turns",    PYONValue([0., 0.5]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_subharmonic_carrier_0_psk_turns",  PYONValue([0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_subharmonic_carrier_1_psk_turns",  PYONValue([0., 0.]), group="{}.psk".format(_argstr))

        # superresolution - custom waveform specification
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
        # convert attenuation from dB to machine units
        self.att_eggs_heating_mu = att_to_mu(self.att_eggs_heating_db * dB)

        # convert freqs to Hz
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # convert build arguments to appropriate values and format as numpy arrays
        freq_eggs_carrier_hz_list = np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_superresolution_sweep_hz_list = np.array(list(self.freq_superresolution_sweep_khz_list)) * kHz
        self.phase_superresolution_sweep_turns_list = np.array(list(self.phase_superresolution_sweep_turns_list))
        self.phase_eggs_heating_ch1_turns_list = np.array(list(self.phase_eggs_heating_ch1_turns_list))
        self.freq_superresolution_osc_base_hz_list = np.array(self.freq_superresolution_osc_khz_list) * kHz + self.freq_global_offset_hz

        # implement variable freq sweep target
        # adjust oscillator phases based on user configuration
        if self.target_freq_sweep == "Secular":
            self.freq_update_arr = np.array([-1., 1., 0., 0])
        elif self.target_freq_sweep == "Probe":
            self.freq_update_arr = np.array([1., 1., 0., 0])
        elif self.target_freq_sweep == "X1":
            self.freq_update_arr = np.array([0., 0., 1., 1.])
        elif self.target_freq_sweep == "X2":
            self.freq_update_arr = np.array([0., 0., -1., 1.])

        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_sweep_turns = np.arange(len(self.phase_superresolution_sweep_turns_list))

        # create config data structure
        self.config_experiment_list = np.zeros((
            len(freq_sideband_readout_ftw_list) *
            len(freq_eggs_carrier_hz_list) *
            len(self.freq_superresolution_sweep_hz_list) *
            len(self.phase_superresolution_sweep_turns_list) *
            len(self.phase_eggs_heating_ch1_turns_list) *
            len(time_readout_mu_list),
        6), dtype=float)
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_experiment_list[:, [1, 2, -3, -2, -1, 0]] = np.stack(np.meshgrid(
            freq_eggs_carrier_hz_list,
            self.freq_superresolution_sweep_hz_list,
            self.waveform_index_to_phase_sweep_turns,
            self.phase_eggs_heating_ch1_turns_list,
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
        if type(self.ampl_superresolution_osc_frac_list) is list:
            if len(self.ampl_superresolution_osc_frac_list) != 4:
                raise ValueError("Error: phaser oscillator amplitude array must have length 4.")
            elif np.sum(self.ampl_superresolution_osc_frac_list) >= 100.:
                raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude array must be a list.")

        if type(self.phase_superresolution_osc_turns_list) is list:
            if len(self.phase_superresolution_osc_turns_list) != 4:
                raise ValueError("Error: phaser oscillator phase array must have length 4.")
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        if type(self.freq_superresolution_osc_khz_list) is not list:
            raise ValueError("Error: phaser oscillator frequency array must be a list.")
        elif len(self.freq_superresolution_osc_khz_list) != 4:
            raise ValueError("Error: phaser oscillator frequency array must have length 4.")
        max_osc_freq_hz = (
                np.max(list(self.freq_superresolution_sweep_khz_list)) * kHz +
                max(self.freq_superresolution_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        min_osc_freq_hz = (
                np.min(list(self.freq_superresolution_sweep_khz_list)) * kHz +
                max(self.freq_superresolution_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        # check that PSK schedule is valid
        num_psk_blocks = len(self.phase_superresolution_rsb_psk_turns)
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (type(psk_schedule) is not list) or (len(psk_schedule) != num_psk_blocks)
            for psk_schedule in (
                self.phase_superresolution_rsb_psk_turns, self.phase_superresolution_bsb_psk_turns,
                self.phase_subharmonic_carrier_0_psk_turns, self.phase_subharmonic_carrier_1_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule. All PSK schedules must be of same length.")

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.waveform_index_to_pulseshaper_vals =   list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_id =     np.zeros(len(self.phase_superresolution_sweep_turns_list),
                                                             dtype=np.int32)   # store pulseshaper waveform ID

        # set up blocks & timing for phaser waveform
        num_blocks = 1
        if self.enable_phase_shift_keying:  num_blocks = len(self.phase_osc0_psk_turns)
        time_block_us = self.time_eggs_heating_us / num_blocks

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = np.zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = np.array(self.ampl_superresolution_osc_frac_list)

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        phase_osc_update_delay_turns_list = (
                (self.freq_superresolution_osc_base_hz_list +
                 self.freq_update_arr * np.mean(self.freq_superresolution_sweep_hz_list)) *
                np.array([0, 40e-9, 80e-9, 80e-9]) # account for relative delay period between each osc
        )
        _osc_vals_blocks[:, :, 1] += np.array(self.phase_superresolution_osc_turns_list) + phase_osc_update_delay_turns_list

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            _osc_vals_blocks[:, :, 1] += np.array([
                self.phase_superresolution_rsb_psk_turns, self.phase_superresolution_bsb_psk_turns,
                self.phase_subharmonic_carrier_0_psk_turns, self.phase_subharmonic_carrier_1_psk_turns
            ][:self._num_phaser_oscs]).transpose()

        # adjust oscillator phases based on user configuration
        if self.target_phase_sweep == "RSB":
            phas_update_arr = np.array([1., 0., 0., 0.])[:self._num_phaser_oscs]
        elif self.target_phase_sweep == "BSB":
            phas_update_arr = np.array([0., 1., 0., 0.])[:self._num_phaser_oscs]
        elif self.target_phase_sweep == "RSB - BSB":
            phas_update_arr = np.array([1., -1., 0., 0.])[:self._num_phaser_oscs]
        elif self.target_phase_sweep == "Carrier 0":
            phas_update_arr = np.array([0., 0., 1., 0.])[:self._num_phaser_oscs]
        elif self.target_phase_sweep == "Carrier 1":
            phas_update_arr = np.array([0., 0., 0., 1.])[:self._num_phaser_oscs]
        elif self.target_phase_sweep == "Carrier 0 - Carrier 1":
            phas_update_arr = np.array([0., 0., 1., -1.])[:self._num_phaser_oscs]

        # record phaser waveforms - one for each phase
        for phase in self.phase_superresolution_sweep_turns_list:
            # create local copy of _osc_vals_blocks and update with target phase
            # note: no need to deep copy b/c it's filled w/immutables
            _osc_vals_blocks_local = np.copy(_osc_vals_blocks)
            _osc_vals_blocks_local[:, :, 1] += phas_update_arr * phase

            # specify sequence as a list of blocks, where each block is a dict
            # note: have to instantiate locally each loop b/c dicts aren't deep copied
            # todo: no pulse shape blocks - maybe? maybe we want now
            _sequence_blocks_local = [
                {
                    "oscillator_parameters": _osc_vals_blocks_local[_idx_block],
                    "config": {
                        "time_us": time_block_us,
                        "pulse_shaping": self.enable_pulse_shaping and ((_idx_block == 0)
                                                                        or (_idx_block == num_blocks - 1)),
                        # note: we DON'T pulse shape each block here
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping and (_idx_block == 0),
                            "pulse_shape_falling": self.enable_pulse_shaping and (_idx_block == num_blocks - 1),
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
        delay_mu(500000)

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
                phase_ch1_turns =   config_vals[4]
                time_readout_mu =   np.int64(config_vals[5])

                # get corresponding phase and waveform ID from the index
                phase_sweep_turns = self.phase_superresolution_sweep_turns_list[phase_sweep_idx]
                waveform_id = self.waveform_index_to_pulseshaper_id[phase_sweep_idx]
                self.core.break_realtime()

                # create frequency update list for oscillators and set phaser frequencies
                freq_update_list = self.freq_superresolution_osc_base_hz_list + freq_sweep_hz * self.freq_update_arr
                self.phaser_eggs.frequency_configure(
                    # carrier frequency (via DUC)
                    carrier_freq_hz - self.freq_global_offset_hz,
                    # oscillator frequencies
                    [freq_update_list[0], freq_update_list[1],
                     freq_update_list[2], freq_update_list[3], 0.],
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

                '''QVSA PULSE'''
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
        # record phaser sequences onto DMA for each RSB phase
        for i in range(len(self.phase_superresolution_sweep_turns_list)):
            # get waveform for given sweep phase
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self.waveform_index_to_pulseshaper_vals[i]

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            self.waveform_index_to_pulseshaper_id[i] = self.pulse_shaper.waveform_record(
                _wav_data_ampl, _wav_data_phas, _wav_data_time
            )
            self.core.break_realtime()

