import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)

from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper


class EGGSHeatingRamsey(LAXExperiment, Experiment):
    """
    Experiment: EGGS Heating Ramsey

    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones in a Ramsey sequence,
    and measure ion temperature via sideband thermometry.
    """
    name = 'EGGS Heating Ramsey'
    kernel_invariants = {
        'config_eggs_heating_list', 'freq_sideband_readout_ftw_list', 'time_readout_mu_list',
        'time_rf_servo_holdoff_mu', 'freq_eggs_carrier_hz_list', 'freq_eggs_secular_hz_list',
        'phase_ramsey_anti_turns_list', 'phase_eggs_heating_ch1_turns_list', 'waveform_index_to_phase_ramsey_turns',
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                                NumberValue(default=100, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",                           BooleanValue(default=True))
        self.setattr_argument("sub_repetitions",                            NumberValue(default=1, precision=0, step=1, min=1, max=500))

        # get subsequences
        self.initialize_subsequence =                                       InitializeQubit(self)
        self.sidebandcool_subsequence =                                     SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =                                  SidebandReadout(self)
        self.readout_subsequence =                                          Readout(self)
        self.rescue_subsequence =                                           RescueIon(self)

        # EGGS RF
        self.setattr_argument("freq_eggs_heating_carrier_mhz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([83.2028]),
                                                                                    CenterScan(83.20175, 0.05, 0.0005, randomize=True),
                                                                                ],
                                                                                global_min=0.005, global_max=4800, global_step=1,
                                                                                unit="MHz", scale=1, precision=6
                                                                            ), group='EGGS_Heating.frequencies')
        self.setattr_argument("freq_eggs_heating_secular_khz_list",         Scannable(
                                                                                default=[
                                                                                    ExplicitScan([777.5]),
                                                                                    CenterScan(777.5, 4, 0.5, randomize=True),
                                                                                    ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                                ],
                                                                                global_min=0, global_max=10000, global_step=1,
                                                                                unit="kHz", scale=1, precision=3
                                                                            ), group='EGGS_Heating.frequencies')

        # EGGS RF - waveform - timing & phase
        self.setattr_argument("time_readout_us_list",                       Scannable(
                                                                                default=[
                                                                                    ExplicitScan([98.8]),
                                                                                    RangeScan(0, 1500, 100, randomize=True),
                                                                                ],
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, precision=5
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("time_eggs_heating_us",                       NumberValue(default=100, precision=2, step=500, min=0.04, max=100000000), group='EGGS_Heating.ramsey')
        self.setattr_argument("phase_eggs_heating_ch1_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 21, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, precision=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_rsb_turns",               NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')
        self.setattr_argument("phase_eggs_heating_bsb_turns",               NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='EGGS_Heating.waveform.time_phase')

        # EGGS RF - waveform - amplitude - general
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=5., precision=1, step=0.5, min=0, max=31.5), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_rsb_pct",                  NumberValue(default=40., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_bsb_pct",                  NumberValue(default=40., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')
        self.setattr_argument("ampl_eggs_heating_carrier_pct",              NumberValue(default=10., precision=2, step=10, min=0.0, max=99), group='EGGS_Heating.waveform.ampl')

        # EGGS RF - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",                       BooleanValue(default=True), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("type_pulse_shape",                           EnumerationValue(['sine_squared', 'error_function'], default='sine_squared'), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",                NumberValue(default=100, precision=1, step=100, min=10, max=100000), group='EGGS_Heating.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",                NumberValue(default=500, precision=0, step=100, min=100, max=2000), group='EGGS_Heating.pulse_shaping')

        # EGGS RF - Ramsey
        self.setattr_argument("enable_ramsey_delay",                        BooleanValue(default=True), group='EGGS_Heating.ramsey')
        self.setattr_argument("time_ramsey_delay_us",                       NumberValue(default=60, precision=2, step=500, min=0.04, max=100000000), group='EGGS_Heating.ramsey')
        self.setattr_argument("target_ramsey_phase",                        EnumerationValue(['RSB', 'BSB', 'Carrier', 'RSB+BSB'], default='Carrier'), group='EGGS_Heating.ramsey')
        self.setattr_argument("phase_ramsey_anti_turns_list",               Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0., 0.5]),
                                                                                    RangeScan(0, 1.0, 11, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, precision=3
                                                                            ), group='EGGS_Heating.ramsey')

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('ttl10')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)
        self.pulse_shaper = PhaserPulseShaper(self, np.array([0., 0., 0.5, 0., 0.]))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SUBSEQUENCE PARAMETERS'''
        # get readout values
        self.freq_sideband_readout_ftw_list =   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.time_readout_mu_list =             np.array([self.core.seconds_to_mu(time_us * us)
                                                          for time_us in self.time_readout_us_list])

        '''EGGS HEATING - TIMING'''
        # add delay time after EGGS pulse to allow RF servo to re-lock
        self.time_rf_servo_holdoff_mu = self.get_parameter("time_rf_servo_holdoff_us", group="eggs",
                                                           conversion_function=us_to_mu)

        '''EGGS HEATING - CONFIG'''
        # convert build arguments to appropriate values and format as numpy arrays
        self.freq_eggs_carrier_hz_list =                        np.array(list(self.freq_eggs_heating_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =                        np.array(list(self.freq_eggs_heating_secular_khz_list)) * kHz
        self.phase_ramsey_anti_turns_list =                     np.array(list(self.phase_ramsey_anti_turns_list))
        self.phase_eggs_heating_ch1_turns_list =                np.array(list(self.phase_eggs_heating_ch1_turns_list))

        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_ramsey_turns =             np.arange(len(self.phase_ramsey_anti_turns_list))

        # create config data structure
        self.config_eggs_heating_list =                         np.zeros((len(self.freq_sideband_readout_ftw_list) *
                                                                          len(self.freq_eggs_carrier_hz_list) *
                                                                          len(self.freq_eggs_secular_hz_list) *
                                                                          len(self.phase_ramsey_anti_turns_list) *
                                                                          len(self.phase_eggs_heating_ch1_turns_list) *
                                                                          len(self.time_readout_mu_list),
                                                                          6), dtype=float)
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_eggs_heating_list[:, [1, 2, -3, -2, -1, 0]] =   np.stack(np.meshgrid(self.freq_eggs_carrier_hz_list,
                                                                                     self.freq_eggs_secular_hz_list,
                                                                                     self.waveform_index_to_phase_ramsey_turns,
                                                                                     self.phase_eggs_heating_ch1_turns_list,
                                                                                     self.time_readout_mu_list,
                                                                                     self.freq_sideband_readout_ftw_list),
                                                                         -1).reshape(-1, 6)

        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config: np.random.shuffle(self.config_eggs_heating_list)
        # precalculate length of configuration list here to reduce run-time overhead
        self.num_configs = len(self.config_eggs_heating_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.waveform_index_to_pulseshaper_vals =   list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_id =     np.zeros(len(self.phase_ramsey_anti_turns_list), dtype=np.int32)   # store pulseshaper waveform ID

        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        self.pulse_shaper._phase_offsets_turns =    np.array([0., 0., 0.5, 0., 0.])

        # set up the spin echo wizard generally
        # note: time_pulse_us divided by num_blocks to split it equally
        self.spinecho_wizard.time_pulse_us =                self.time_eggs_heating_us
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           True
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        self.enable_ramsey_delay
        self.spinecho_wizard.time_delay_spinecho_us =       self.time_ramsey_delay_us

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence
        _sequence_blocks = np.zeros((2, 3, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_eggs_heating_rsb_pct
        _sequence_blocks[:, 1, 0] = self.ampl_eggs_heating_bsb_pct
        _sequence_blocks[:, 2, 0] = self.ampl_eggs_heating_carrier_pct

        # set rsb & bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        phase_bsb_update_delay_turns = np.mean(self.freq_eggs_secular_hz_list) * (self.phaser_eggs.t_sample_mu * ns)
        _sequence_blocks[:, 0, 1] = self.phase_eggs_heating_rsb_turns
        _sequence_blocks[:, 1, 1] = self.phase_eggs_heating_bsb_turns + phase_bsb_update_delay_turns

        # get ramsey phase target so we can Ramsey on different oscillators
        ramsey_osc_target = 0
        if self.target_ramsey_phase == 'RSB':
            ramsey_osc_target = 0
        elif self.target_ramsey_phase == 'BSB':
            ramsey_osc_target = 1
        elif self.target_ramsey_phase == 'Carrier':
            ramsey_osc_target = 2


        # record EGGS pulse waveforms
        for i in range(len(self.phase_ramsey_anti_turns_list)):
            # update sequence block with ramsey phase
            phase_ramsey_turns = self.phase_ramsey_anti_turns_list[i]
            _sequence_blocks[1, ramsey_osc_target, 1] = phase_ramsey_turns

            # set phase shift for RSB+BSB case
            if self.target_ramsey_phase == 'RSB+BSB':
                _sequence_blocks[1, 0, 1] = phase_ramsey_turns
                _sequence_blocks[1, 1, 1] = phase_ramsey_turns

            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_blocks
            self.spinecho_wizard.calculate_pulseshape()
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_index_to_pulseshaper_vals.append(self.spinecho_wizard.get_waveform())

        # tmp remove
        # _wav_print_idk = self.waveform_index_to_pulseshaper_vals[0]
        # print(_wav_print_idk[0])
        # print(_wav_print_idk[1])
        # print(_wav_print_idk[2])
        # print(_sequence_blocks)
        # self.spinecho_wizard.display_waveform()
        # tmp remove

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions * len(self.config_eggs_heating_list),
                7)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        ### PHASER INITIALIZATION ###
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(31.5 * dB)

        # ensure INT HOLD on RF servo is deactivated
        self.ttl10.off()

        # reset debug triggers
        self.ttl8.off()
        self.ttl9.off()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        # load waveform DMA handles
        self.pulse_shaper.waveform_load()
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # tmp remove
        # set phaser attenuators here to prevent turn on glitch (problematic when cables direct to trap)
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)
        self.core.break_realtime()
        # tmp remove

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # implement sub-repetitions here to avoid initial overhead
            _subrep_iter = 0
            _config_iter = 0

            # sweep experiment configurations
            while _config_iter < self.num_configs:

                '''CONFIGURE'''
                config_vals = self.config_eggs_heating_list[_config_iter]
                # extract values from config list
                freq_readout_ftw =      np.int32(config_vals[0])
                carrier_freq_hz =       config_vals[1]
                sideband_freq_hz =      config_vals[2]
                phase_ramsey_index =    np.int32(config_vals[3])
                phase_ch1_turns =       config_vals[4]
                time_readout_mu =       np.int64(config_vals[5])

                # get corresponding RSB phase and waveform ID from the index
                phase_ramsey_turns = self.phase_ramsey_anti_turns_list[phase_ramsey_index]
                waveform_id = self.waveform_index_to_pulseshaper_id[phase_ramsey_index]
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz, phase_ch1_turns)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                '''EGGS HEATING'''
                self.phaser_run(waveform_id)

                '''READOUT'''
                self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        counts,
                        carrier_freq_hz,
                        sideband_freq_hz,
                        phase_ramsey_turns,
                        phase_ch1_turns,
                        time_readout_mu
                    )
                    self.core.break_realtime()

                '''LOOP CLEANUP'''
                # resuscitate ion
                self.rescue_subsequence.resuscitate()

                # death detection
                self.rescue_subsequence.detect_death(counts)
                self.core.break_realtime()

                # check termination more frequently in case reps are low
                if (_loop_iter % 50) == 0:
                    self.check_termination()
                    self.core.break_realtime()
                _loop_iter += 1

                # handle sub-repetition logic
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
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

        '''CLEANUP'''
        self.core.break_realtime()
        self.ttl10.off()


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
        # activate integrator hold
        self.ttl10.on()

        # # set phaser attenuators - warning creates turn on glitch
        # at_mu(self.phaser_eggs.get_next_frame_mu())
        # self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
        # delay_mu(self.phaser_eggs.t_sample_mu)
        # self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

        # tmp remove
        # add delay time after integrator hold to reduce effect of turn-on glitches
        delay_mu(self.time_rf_servo_holdoff_mu)
        # tmp remove

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        with parallel:
            # activate debug TTLs
            self.ttl8.on()

            # play recorded waveform
            self.pulse_shaper.waveform_playback(waveform_id)

        # EGGS - STOP
        # stop all oscillators (doesn't unset attenuators)
        self.phaser_eggs.phaser_stop()
        # deactivate integrator hold
        self.ttl10.off()
        # stop debug TTLs
        self.ttl8.off()
        # add delay time after EGGS pulse to allow RF servo to re-lock
        delay_mu(self.time_rf_servo_holdoff_mu)

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each ramsey phase
        for i in range(len(self.phase_ramsey_anti_turns_list)):

            # get waveform for given ramsey phase
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self.waveform_index_to_pulseshaper_vals[i]
            self.core.break_realtime()

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1 ms)
            _wav_idx = self.pulse_shaper.waveform_record(_wav_data_ampl, _wav_data_phas, _wav_data_time)
            self.waveform_index_to_pulseshaper_id[i] = _wav_idx
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, phase_ch1_offset_turns: TFloat) -> TNone:
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the carrier frequency (in Hz).
            sideband_freq_hz        (float)     : the sideband frequency (in Hz).
            phase_ch1_offset_turns  (float)     : the phase offset for CH1 (in turns).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        phase_ch1_turns = phase_ch1_offset_turns + (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns)

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(phase_ch1_turns)
        self.phaser_eggs.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)


    '''
    ANALYSIS
    '''
    def analyze(self):
        # should be backward-compatible with experiment which had no sub-repetitions
        try:
            sub_reps = self.sub_repetitions
        except Exception as e:
            sub_reps = 1

        # handle errors from data processing
        try:
            # print results
            print("\tResults - EGGS Heating:")

            # convert dataset to array
            dataset = np.array(self.results)

            ## determine scan type
            # carrier sweep
            if len(np.unique(dataset[:, 2])) > 1:       sorting_col_num = 2
            # secular frequency
            elif len(np.unique(dataset[:, 3])) > 1:     sorting_col_num = 3
            else:                                       sorting_col_num = 0

            ## convert results to sideband ratio
            ratios, ave_rsb, ave_bsb, std_rsb, std_bsb, scanning_freq_MHz = extract_ratios(dataset, sorting_col_num,
                                                                                           1, 0,
                                                                                           self.repetitions, sub_reps)
            phonons = convert_ratios_to_coherent_phonons(ratios)

            ## process secular frequency sweep
            if sorting_col_num == 3:
                fit_params_secular, fit_err_secular, _ = fitSincGeneric(scanning_freq_MHz, phonons)
                phonon_n = fit_params_secular[0]
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_secular',  fit_params_secular)
                self.set_dataset('fit_err_secular',     fit_err_secular)

                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_secular, fit_err_secular]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

                # print results to log
                print("\t\tSecular: {:.4f} +/- {:.5f} kHz".format(fit_params_secular[1] * 1e3, fit_err_secular[1] * 1e3))


            ## process sideband readout sweep
            else:
                rsb_freqs_MHz, bsb_freqs_MHz, _ =       extract_sidebands_freqs(scanning_freq_MHz)
                fit_params_rsb, fit_err_rsb, fit_rsb =  fitSincGeneric(rsb_freqs_MHz, ave_rsb)
                fit_params_bsb, fit_err_bsb, fit_bsb =  fitSincGeneric(bsb_freqs_MHz, ave_bsb)
                phonon_n = fit_params_rsb[0] / (fit_params_bsb[0] - fit_params_rsb[0])
                # todo: implement
                phonon_err = 0

                # save results to hdf5 as a dataset
                self.set_dataset('fit_params_rsb',  fit_params_rsb)
                self.set_dataset('fit_params_bsb',  fit_params_bsb)
                self.set_dataset('fit_err_rsb',     fit_err_rsb)
                self.set_dataset('fit_err_bsb',     fit_err_bsb)

                # save results to dataset manager for dynamic experiments
                res_dj = [[phonon_n, phonon_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
                self.set_dataset('temp.eggsheating.results', res_dj, broadcast=True, persist=False, archive=False)
                self.set_dataset('temp.eggsheating.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

                # print results to log
                print("\t\tRSB: {:.4f} +/- {:.5f}".format(float(fit_params_rsb[1]) / 2., float(fit_err_rsb[1]) / 2.))
                print("\t\tBSB: {:.4f} +/- {:.5f}".format(float(fit_params_bsb[1]) / 2., float(fit_err_bsb[1]) / 2.))

        except Exception as e:
            print("Warning: unable to process data.")
            print(repr(e))
