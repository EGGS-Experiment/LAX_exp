import numpy as np
from sipyco import pyon
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM, SidebandReadout, QubitRAP
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper


class EGGSQLS(LAXExperiment, Experiment):
    """
    Experiment: EGGS QLS

    Single-tone EGGS-based Quantum Logic Spectroscopy.
    Drive a motional BSB with a single tone, then apply a carrier pulse to reinitialize the state.
    """
    name = 'EGGS QLS'
    kernel_invariants = {
        # hardware parameters
        'att_eggs_qls_mu', 'freq_eggs_secular_hz_list', 'phase_eggs_qls_ch1_turns_list', 'num_configs',
        'waveform_index_to_pulseshaper_vals',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'enable_RAP',

        # configs
        'profile_729_readout', 'profile_729_SBC', 'profile_729_RAP', 'config_experiment_list'
    }

    def build_experiment(self):
        # exp-specific variables
        _argstr = "QLS"              # string to use for arguments

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=80, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config", BooleanValue(default=False))
        self.setattr_argument("sub_repetitions", NumberValue(default=1, precision=0, step=1, min=1, max=500))
        self.setattr_argument("readout_type", EnumerationValue(["Sideband Ratio", "RAP"], default="RAP"))

        # allocate relevant beam profiles
        self.profile_729_readout = 0
        self.profile_729_SBC = 1
        self.profile_729_RAP = 2

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )
        self.sidebandreadout_subsequence = SidebandReadout(self, profile_dds=self.profile_729_readout)
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # Sideband Readout - extra argument
        self.setattr_argument("time_readout_us_list", Scannable(
            default=[
                ExplicitScan([120.5]),
                RangeScan(0, 1500, 100, randomize=True),
            ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ), group='sideband_readout')

        # RAP-based readout
        self.setattr_argument("att_rap_db", NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit='dB', scale=1.),
                              group="RAP")
        self.setattr_argument("ampl_rap_pct", NumberValue(default=50., precision=3, step=5, min=1, max=50, unit='%', scale=1.),
                              group="RAP")
        self.setattr_argument("freq_rap_center_mhz", NumberValue(default=101.3318, precision=6, step=1e-2, min=60, max=200, unit='MHz', scale=1.),
                              group='RAP')
        self.setattr_argument("freq_rap_dev_khz",   NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit='kHz', scale=1.),
                              group='RAP')
        self.setattr_argument("time_rap_us", NumberValue(default=500., precision=3, min=1, max=1e5, step=1, unit='us', scale=1.),
                              group="RAP")

        # EGGS RF
        self.setattr_argument("freq_eggs_qls_carrier_mhz_list", Scannable(
            default=[
                ExplicitScan([80.]),
                CenterScan(83.20175, 0.05, 0.0005, randomize=False),
                RangeScan(82, 84, 100, randomize=False),
            ],
            global_min=0.005, global_max=4800, global_step=1,
            unit="MHz", scale=1, precision=6
        ), group='{}.frequencies'.format(_argstr))
        self.setattr_argument("freq_eggs_qls_secular_khz_list", Scannable(
            default=[
                ExplicitScan([1303]),
                CenterScan(1303, 4, 0.1, randomize=True),
                # ExplicitScan([767.2, 319.2, 1582, 3182]),
            ],
            global_min=0, global_max=10000, global_step=1,
            unit="kHz", scale=1, precision=3
        ), group='{}.frequencies'.format(_argstr))

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_eggs_qls_us", NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000, unit='us', scale=1.),
                              group='{}.waveform.time_phase'.format(_argstr))
        self.setattr_argument("time_eggs_carrier_us", NumberValue(default=100, precision=2, step=500, min=0.04, max=100000000, unit='us', scale=1.),
                              group='{}.waveform.time_phase'.format(_argstr))

        self.setattr_argument("phase_eggs_qls_ch1_turns_list", Scannable(
            default=[
                ExplicitScan([0.4668]),
                RangeScan(0, 1.0, 2, randomize=True),
            ],
            global_min=-1.0, global_max=1.0, global_step=0.01,
            unit="turns", scale=1, precision=3
        ), group='{}.waveform.time_phase'.format(_argstr))
        self.setattr_argument("phase_eggs_qls_bsb_turns_list", Scannable(
            default=[
                ExplicitScan([0]),
                RangeScan(0, 1.0, 2, randomize=True),
            ],
            global_min=-1.0, global_max=1.0, global_step=0.01,
            unit="turns", scale=1, precision=3
        ), group='{}.waveform.time_phase'.format(_argstr))

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("att_eggs_qls_db", NumberValue(default=10., precision=1, step=0.5, min=0, max=31.5, unit='dB', scale=1.),
                              group='{}.waveform.ampl'.format(_argstr))
        self.setattr_argument("ampl_eggs_qls_bsb_pct", NumberValue(default=40., precision=2, step=10, min=0.0, max=99, unit='%', scale=1.),
                              group='{}.waveform.ampl'.format(_argstr))
        self.setattr_argument("ampl_eggs_qls_carrier_pct", NumberValue(default=40, precision=2, step=10, min=0.0, max=99, unit='%', scale=1.),
                              group='{}.waveform.ampl'.format(_argstr))

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_qls_pulse_shaping", BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_qls_pulse_shape", EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_qls_pulse_shape_rolloff_us", NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit='us', scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_qls_pulse_shape_sample_khz", NumberValue(default=1000, precision=0, step=100, min=100, max=5000, unit='us', scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        self.setattr_argument("enable_carrier_pulse_shaping", BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_carrier_pulse_shape", EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_carrier_pulse_shape_rolloff_us", NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit='us', scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_carrier_pulse_shape_sample_khz", NumberValue(default=1000, precision=0, step=100, min=100, max=5000, unit='kHz', scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        # instantiate RAP here
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )
        # get relevant devices & instantiate helper objects
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')
        self.spinecho_wizard = SpinEchoWizardRDX(self)
        # phase delays set field geometries (0.5 for osc_1 => dipole)
        self.pulse_shaper = PhaserPulseShaper(self, np.array([0., 0.5, 0., 0., 0.]))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE/VALIDATE INPUTS & CHECK ERRORS'''
        self._prepare_argument_checks()

        '''SUBSEQUENCE PARAMETERS'''
        # prepare RAP arguments
        self.att_rap_mu = att_to_mu(self.att_rap_db * dB)
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        # configure readout method
        if self.readout_type == 'RAP':
            self.enable_RAP = True
            freq_sideband_readout_ftw_list = np.array([self.freq_rap_center_ftw], dtype=np.int32)
            time_readout_mu_list = np.array([self.time_rap_mu], dtype=np.int64)
        elif self.readout_type == 'Sideband Ratio':
            self.enable_RAP = False
            freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
            time_readout_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                             for time_us in self.time_readout_us_list])
        else:
            raise ValueError("Invalid readout type. Must be one of (Sideband Ratio, RAP).")

        '''EGGS QLS - CONFIG'''
        # convert attenuation from dB to machine units
        self.att_eggs_qls_mu = att_to_mu(self.att_eggs_qls_db * dB)

        # convert build arguments to appropriate values and format as numpy arrays
        freq_eggs_carrier_hz_list = np.array(list(self.freq_eggs_qls_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list = np.array(list(self.freq_eggs_qls_secular_khz_list)) * kHz
        self.phase_eggs_qls_bsb_turns_list = np.array(list(self.phase_eggs_qls_bsb_turns_list))
        self.phase_eggs_qls_ch1_turns_list = np.array(list(self.phase_eggs_qls_ch1_turns_list))

        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_bsb_turns = np.arange(len(self.phase_eggs_qls_bsb_turns_list))

        # create config data structure
        self.config_experiment_list = np.zeros((
            len(freq_sideband_readout_ftw_list) *
            len(freq_eggs_carrier_hz_list) *
            len(self.freq_eggs_secular_hz_list) *
            len(self.phase_eggs_qls_ch1_turns_list) *
            len(self.phase_eggs_qls_bsb_turns_list) *
            len(time_readout_mu_list),
            6), dtype=float)
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_experiment_list[:, [1, 2, -3, -2, -1, 0]] = np.stack(np.meshgrid(
            freq_eggs_carrier_hz_list,
            self.freq_eggs_secular_hz_list,
            self.waveform_index_to_phase_bsb_turns,
            self.phase_eggs_qls_ch1_turns_list,
            time_readout_mu_list,
            freq_sideband_readout_ftw_list),
            -1).reshape(-1, 6)

        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config:   np.random.shuffle(self.config_experiment_list)
        # precalculate length of configuration list here to reduce run-time overhead
        self.num_configs = len(self.config_experiment_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for pulse waveforms
        num_blocks = 2
        self.waveform_index_to_pulseshaper_vals = list()  # store compiled waveforms
        self.waveform_index_to_pulseshaper_id = np.zeros(len(self.phase_eggs_qls_bsb_turns_list),
                                                         dtype=np.int32)  # store pulseshaper waveform ID

        '''DESIGN WAVEFORM SEQUENCE'''
        # create separate bare waveform block sequences for CH0 and CH1
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        _osc_vals_blocks = np.zeros((num_blocks, 2, 2), dtype=float)
        _osc_vals_blocks[0, 0, 0] = self.ampl_eggs_qls_bsb_pct
        _osc_vals_blocks[1, 1, 0] = self.ampl_eggs_qls_carrier_pct

        # note: no phase delay needed for carrier or bsb

        '''COMPILE WAVEFORM SEQUENCE'''
        for i, phase in enumerate(self.phase_eggs_qls_bsb_turns_list):
            # create local copy of _sequence_blocks
            # note: no need to deep copy b/c it's filled w/immutables
            # note: have to obtain different copies so they don't point to same object and overwrite it
            _osc_vals_blocks_local = np.copy(_osc_vals_blocks)
            _osc_vals_blocks_local[0, 0, 1] += phase

            _sequence_blocks_local = [
                {
                    "oscillator_parameters": _osc_vals_blocks_local[0, :, :],
                    "config": {
                        "time_us": self.time_eggs_qls_us,
                        "pulse_shaping": self.enable_qls_pulse_shaping,
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_qls_pulse_shape,
                            "pulse_shape_rising": self.enable_qls_pulse_shaping,
                            "pulse_shape_falling": self.enable_qls_pulse_shaping,
                            "sample_rate_khz": self.freq_qls_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_qls_pulse_shape_rolloff_us
                        }
                    }
                },
                {
                    "oscillator_parameters": _osc_vals_blocks_local[1, :, :],
                    "config": {
                        "time_us": self.time_eggs_carrier_us,
                        "pulse_shaping": self.enable_carrier_pulse_shaping,
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_carrier_pulse_shape,
                            "pulse_shape_rising": self.enable_carrier_pulse_shaping,
                            "pulse_shape_falling": self.enable_carrier_pulse_shaping,
                            "sample_rate_khz": self.freq_carrier_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_carrier_pulse_shape_rolloff_us
                        }
                    }
                }
            ]

            # compile waveform and store in holding structure
            self.waveform_index_to_pulseshaper_vals.append(self.spinecho_wizard.compile_waveform(_sequence_blocks_local))

    def _prepare_argument_checks(self):
        """
        Check experiment arguments for validity.
        """
        # ensure phaser amplitudes sum to less than 100%
        total_phaser_channel_amplitude = (self.ampl_eggs_qls_bsb_pct +
                                          self.ampl_eggs_qls_carrier_pct)
        if total_phaser_channel_amplitude > 100.:
            raise ValueError("Error: phaser oscillator amplitudes exceed 100%.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_eggs_qls_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions * len(self.config_experiment_list),
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
        _loop_iter = 0 # used to check_termination more frequently

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            # implement sub-repetitions here to avoid initial overhead
            _subrep_iter = 0
            _config_iter = 0

            # sweep experiment configurations
            while _config_iter < self.num_configs:

                '''CONFIGURE'''
                # extract values from config list
                config_vals = self.config_experiment_list[_config_iter]
                freq_readout_ftw =  np.int32(config_vals[0])
                carrier_freq_hz =   config_vals[1]
                sideband_freq_hz =  config_vals[2]
                phase_bsb_index =   np.int32(config_vals[3])
                phase_ch1_turns =   config_vals[4]
                time_readout_mu =   np.int64(config_vals[5])
                phase_bsb_turns =   self.phase_eggs_qls_bsb_turns_list[phase_bsb_index]
                waveform_id =       self.waveform_index_to_pulseshaper_id[phase_bsb_index]

                # configure EGGS tones and set readout frequency
                self.core.break_realtime()
                self.phaser_eggs.frequency_configure(carrier_freq_hz,
                                                     [sideband_freq_hz, 0., 0., 0., 0.],
                                                     phase_ch1_turns)
                delay_mu(25000)
                if not self.enable_RAP:
                    self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                      profile=self.profile_729_readout, phase_mode=PHASE_MODE_CONTINUOUS)
                    delay_mu(25000)

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state & sideband cool to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                '''EGGS HEATING'''
                self.phaser_run(waveform_id)

                '''READOUT'''
                # readout via RAP or sideband ratio
                if self.enable_RAP:
                    self.qubit.set_att_mu(self.att_rap_mu)
                    self.rap_subsequence.run_rap(time_readout_mu)
                else:
                    self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()

                # get counts & clean up shot
                counts = self.readout_subsequence.fetch_count()
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(counts)

                '''LOOP CLEANUP'''
                # update dataset
                self.update_results(
                    freq_readout_ftw,
                    counts,
                    carrier_freq_hz,
                    sideband_freq_hz,
                    phase_bsb_turns,
                    phase_ch1_turns,
                    time_readout_mu
                )

                # check termination more frequently in case reps are low
                if _loop_iter % 50 == 0:
                    self.check_termination()
                _loop_iter += 1

                # handle sub-repetition logic
                # handle sub-repetition logic
                # todo: document
                # todo: fix - when we do RAP, we don't have to do subrep stuff
                if _config_iter % 2 == 1:
                    _subrep_iter += 1
                    if _subrep_iter < self.sub_repetitions:
                        _config_iter -= 1
                    else:
                        _subrep_iter = 0
                        _config_iter += 1
                else:
                    _config_iter += 1

            # rescue ion as needed
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
        self.phaser_eggs.phaser_setup(self.att_eggs_qls_mu, self.att_eggs_qls_mu)

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
        for i in range(len(self.phase_eggs_qls_bsb_turns_list)):
            # get waveform for given parameters
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self._get_compiled_waveform(i)

            # record phaser pulse sequence and save returned waveform ID
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
        """
        return self.waveform_index_to_pulseshaper_vals[wav_idx]


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
        ch1_sweep_bool = len(self.phase_eggs_qls_ch1_turns_list) > 1
        turns_sweep_bool = ch1_sweep_bool

        # print results
        print("\tResults - EGGS Heating:")

        # sweep over ch1_turns
        for ch1_turns in self.phase_eggs_qls_ch1_turns_list:

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
                    ccb_command = f'$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_qls_RDX_{ch1_turns}'
                    group = 'plotting.eggs_qls.ch1_sweep'
                    dataset_name = f'temp.plotting.results_eggs_qls_RDX_{ch1_turns}'
                    applet_name = f"EGGS Heating - RDX - CH1 Turns: {ch1_turns}"
                else:
                    ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_qls_RDX'
                    group = 'plotting.eggs_qls'
                    dataset_name = f'temp.plotting.results_eggs_qls_RDX'
                    applet_name = f"EGGS Heating - RDX"

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

                    # determine group to plot applet in
                    if sorting_col_num == 3:
                        group += '.secular'
                    else:
                        group += '.carrier'

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

                    group += '.sidebands'

                else:
                    # if unknown quantity was scanned ignore analysis
                    plotting_results = {}
                    res_dj = None
                    group = None

                # record values in dataset
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False,
                                 archive=False)
                self.set_dataset(dataset_name, pyon.encode(plotting_results), broadcast=True)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False,
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

                self.set_dataset('temp.plotting.results_eggs_qls_RDX_ch1_sweep',
                                 pyon.encode(plotting_results), broadcast=True)

                self.ccb.issue("create_applet", f"EGGS Heating - RDX",
                               '$python -m LAX_exp.applets.plot_matplotlib '
                               'temp.plotting.results_eggs_qls_RDX_ch1_sweep'
                               ' --num-subplots 1', group=['plotting', 'eggs_qls', 'ch1_sweep'])


            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))
