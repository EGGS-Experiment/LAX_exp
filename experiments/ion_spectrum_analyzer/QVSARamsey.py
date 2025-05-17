import numpy as np
from sipyco import pyon
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuousRAM, SidebandReadout)

from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper


class QVSARamsey(LAXExperiment, Experiment):
    """
    Experiment: QVSA Ramsey

    Cool the ions to the ground motional state via sideband cooling,
    then apply a motional Ramsey sequence using the QVSA effect,
    then measure the residual ion motion using either sideband thermometry or rabi divination.
    """
    name = 'QVSA Ramsey'
    kernel_invariants = {
        # hardware parameters
        'att_qvsa_mu', 'freq_eggs_secular_hz_list',
        'phase_ramsey_anti_turns_list', 'phase_qvsa_ch1_turns_list', 'waveform_index_to_phase_ramsey_turns',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence',
        'rescue_subsequence',
        
        # configs
        'profile_729_readout', 'profile_729_SBC', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",        NumberValue(default=100, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("enable_linetrigger", BooleanValue(default=True))

        # allocate relevant beam profiles
        self.profile_729_readout =  0
        self.profile_729_SBC =      1

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=500
        )
        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_readout)
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # Sideband Readout - extra argument
        self.setattr_argument("time_readout_us_list",   Scannable(
                                                                    default=[
                                                                        ExplicitScan([50.8]),
                                                                        RangeScan(0, 1500, 100, randomize=True),
                                                                    ],
                                                                    global_min=1, global_max=100000, global_step=1,
                                                                    unit="us", scale=1, precision=5
                                                                ), group='sideband_readout')

        # QVSA - frequency
        self.setattr_argument("freq_qvsa_carrier_mhz_list", Scannable(
                                                                            default=[
                                                                                ExplicitScan([80.2028]),
                                                                                CenterScan(83.20175, 0.05, 0.0005, randomize=True),
                                                                            ],
                                                                            global_min=0.005, global_max=4800, global_step=1,
                                                                            unit="MHz", scale=1, precision=6
                                                                        ), group='QVSA.frequencies')
        self.setattr_argument("freq_qvsa_secular_khz_list", Scannable(
                                                                            default=[
                                                                                ExplicitScan([777.5]),
                                                                                CenterScan(777.5, 4, 0.5, randomize=True),
                                                                                ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                            ],
                                                                            global_min=0, global_max=10000, global_step=1,
                                                                            unit="kHz", scale=1, precision=3
                                                                        ), group='QVSA.frequencies')

        # QVSA - Ramsey
        self.setattr_argument("enable_ramsey_delay",    BooleanValue(default=True), group='QVSA.ramsey')
        self.setattr_argument("time_ramsey_delay_us",   NumberValue(default=60, precision=2, step=500, min=0.04, max=100000000), group='QVSA.ramsey')
        self.setattr_argument("target_ramsey_phase",    EnumerationValue(['RSB', 'BSB', 'Carrier', 'RSB+BSB'], default='RSB+BSB'), group='QVSA.ramsey')
        self.setattr_argument("phase_ramsey_anti_turns_list",   Scannable(
                                                                        default=[
                                                                            ExplicitScan([0., 0.5]),
                                                                            RangeScan(0, 1.0, 11, randomize=True),
                                                                        ],
                                                                        global_min=0.0, global_max=1.0, global_step=1,
                                                                        unit="turns", scale=1, precision=3
                                                                    ), group='QVSA.ramsey')

        # QVSA - waveform
        self.setattr_argument("att_qvsa_db",                NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5), group='QVSA.waveform')
        self.setattr_argument("time_qvsa_us",               NumberValue(default=100, precision=2, step=500, min=0.04, max=100000000), group='QVSA.waveform')
        self.setattr_argument("ampl_qvsa_osc_frac_list",    PYONValue([10., 10., 0.]), group="QVSA.waveform")
        self.setattr_argument("phase_qvsa_osc_turns_list",  PYONValue([0., 0., 0.]), group="QVSA.waveform")
        self.setattr_argument("phase_qvsa_ch1_turns_list",  Scannable(
                                                                        default=[
                                                                            ExplicitScan([0.2]),
                                                                            RangeScan(0, 1.0, 21, randomize=True),
                                                                        ],
                                                                        global_min=0.0, global_max=1.0, global_step=1,
                                                                        unit="turns", scale=1, precision=3
                                                                    ), group='QVSA.waveform')

        # QVSA - waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=True), group="QVSA.waveform.pulse_shaping")
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group="QVSA.waveform.pulse_shaping")
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=50, precision=1, step=100, min=0.2, max=100000),
                              group="QVSA.waveform.pulse_shaping")
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1000, precision=0, step=100, min=100, max=5000),
                              group="QVSA.waveform.pulse_shaping")

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')
        self.setattr_device('trigger_line')
        # tmp remove
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        # tmp remove

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)
        self.pulse_shaper = PhaserPulseShaper(self, np.array([0., 0., 0.5, 0., 0.]))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE & VALIDATE INPUTS'''
        self._prepare_argument_checks()

        '''SUBSEQUENCE PARAMETERS'''
        # get readout values
        freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        time_readout_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                         for time_us in self.time_readout_us_list])

        '''QVSA - CONFIG'''
        # convert attenuation from dB to machine units
        self.att_qvsa_mu = att_to_mu(self.att_qvsa_db * dB)

        # convert build arguments to appropriate values and format as numpy arrays
        freq_qvsa_carrier_hz_list =             np.array(list(self.freq_qvsa_carrier_mhz_list)) * MHz
        self.freq_eggs_secular_hz_list =        np.array(list(self.freq_qvsa_secular_khz_list)) * kHz
        self.phase_ramsey_anti_turns_list =     np.array(list(self.phase_ramsey_anti_turns_list))
        self.phase_qvsa_ch1_turns_list =        np.array(list(self.phase_qvsa_ch1_turns_list))

        # map phase to index to facilitate waveform recording
        self.waveform_index_to_phase_ramsey_turns = np.arange(len(self.phase_ramsey_anti_turns_list))

        # create config data structure
        self.config_experiment_list = np.zeros(
            (len(freq_sideband_readout_ftw_list) * len(freq_qvsa_carrier_hz_list) *
             len(self.freq_eggs_secular_hz_list) * len(self.phase_ramsey_anti_turns_list) *
             len(self.phase_qvsa_ch1_turns_list) * len(time_readout_mu_list),
             6),
            dtype=float
        )
        # note: sideband readout frequencies are at the end of the meshgrid
        # to ensure successive rsb/bsb measurements are adjacent
        self.config_experiment_list[:, [1, 2, -3, -2, -1, 0]] = np.stack(
            np.meshgrid(freq_qvsa_carrier_hz_list,
                        self.freq_eggs_secular_hz_list,
                        self.waveform_index_to_phase_ramsey_turns,
                        self.phase_qvsa_ch1_turns_list,
                        time_readout_mu_list,
                        freq_sideband_readout_ftw_list),
            -1).reshape(-1, 6)
        # always randomize config
        np.random.shuffle(self.config_experiment_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # check that input amplitude/phase arrays are valid
        if type(self.ampl_qvsa_osc_frac_list) is list:
            if len(self.ampl_qvsa_osc_frac_list) != 3:
                raise ValueError("Error: phaser oscillator amplitude array must have length 4.")
            elif np.sum(self.ampl_qvsa_osc_frac_list) >= 100.:
                raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude array must be a list.")

        if type(self.phase_qvsa_osc_turns_list) is list:
            if len(self.phase_qvsa_osc_turns_list) != 3:
                raise ValueError("Error: phaser oscillator phase array must have length 4.")
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        max_osc_freq_hz = max(list(self.freq_qvsa_secular_khz_list)) * kHz
        min_osc_freq_hz = min(list(self.freq_qvsa_secular_khz_list)) * kHz
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_output_freqs_hz = np.array(list(self.freq_qvsa_carrier_mhz_list)) * MHz
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - min(phaser_output_freqs_hz))
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - max(phaser_output_freqs_hz))
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA pulse.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for QVSA pulse waveforms
        self.waveform_index_to_pulseshaper_vals =   list()      # store compiled waveforms
        self.waveform_index_to_pulseshaper_id =     np.zeros(len(self.phase_ramsey_anti_turns_list), dtype=np.int32)   # store pulseshaper waveform ID

        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        self.pulse_shaper._phase_offsets_turns =    np.array([0., 0., 0.5, 0., 0.])

        # set up the spin echo wizard generally
        # note: time_pulse_us is amount of time for each block
        self.spinecho_wizard.time_pulse_us =                self.time_qvsa_us
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           True
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        self.enable_ramsey_delay
        self.spinecho_wizard.time_delay_spinecho_us =       self.time_ramsey_delay_us

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence & set amplitudes
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        _sequence_blocks = np.zeros((2, 3, 2), dtype=float)
        _sequence_blocks[:, :, 0] = np.array(self.ampl_qvsa_osc_frac_list)

        # set rsb & bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        phase_bsb_update_delay_turns = np.mean(self.freq_eggs_secular_hz_list) * (self.phaser_eggs.t_sample_mu * ns)
        _sequence_blocks[:, :, 1] += np.array(self.phase_qvsa_osc_turns_list)
        _sequence_blocks[:, 1, 1] += phase_bsb_update_delay_turns

        # get ramsey phase target so we can Ramsey on different oscillators
        if self.target_ramsey_phase == 'RSB':
            phas_update_arr = np.array([1., 0., 0.])
        elif self.target_ramsey_phase == 'BSB':
            phas_update_arr = np.array([0., 1., 0.])
        elif self.target_ramsey_phase == 'Carrier':
            phas_update_arr = np.array([0., 0., 1.])
        elif self.target_ramsey_phase == 'RSB+BSB':
            phas_update_arr = np.array([1., 1., 0.])

        # record QVSA pulse waveforms
        for i, phase in enumerate(self.phase_ramsey_anti_turns_list):
            # create local copy of _sequence_blocks and update with ramsey phase
            # note: no need to deep copy b/c it's filled w/immutables
            _sequence_blocks_local = np.copy(_sequence_blocks)
            _sequence_blocks_local[1, :, 1] += phas_update_arr * phase

            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_blocks_local
            self.spinecho_wizard.calculate_pulseshape()
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_index_to_pulseshaper_vals.append(self.spinecho_wizard.get_waveform())

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

            # sweep experiment configurations
            for config_vals in self.config_experiment_list:

                # tmp remove
                # turn on rescue beams while waiting
                self.pump.rescue()
                self.repump_cooling.on()
                self.repump_qubit.on()
                self.pump.on()
                # tmp remove

                '''CONFIGURE'''
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
                self.phaser_eggs.frequency_configure(carrier_freq_hz,
                                                     [-sideband_freq_hz, sideband_freq_hz, 0., 0., 0.],
                                                     phase_ch1_turns)
                delay_mu(10000)
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_readout)
                delay_mu(10000)

                # wait for linetrigger
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, self.trigger_line.time_holdoff_mu)

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                '''QVSA RAMSEY PULSE'''
                self.phaser_run(waveform_id)

                '''READOUT'''
                self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
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
                # resuscitate ion & detect death
                self.rescue_subsequence.resuscitate()
                self.rescue_subsequence.detect_death(counts)
                self.core.break_realtime()

                # check termination more frequently in case reps are low
                if (_loop_iter % 50) == 0:
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
        Run the main QVSA pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        """
        # QVSA - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_qvsa_mu, self.att_qvsa_mu)

        # QVSA - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # QVSA - STOP
        # stop all output & clean up hardware (e.g. amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

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
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            _wav_idx = self.pulse_shaper.waveform_record(_wav_data_ampl, _wav_data_phas, _wav_data_time)
            self.waveform_index_to_pulseshaper_id[i] = _wav_idx
            self.core.break_realtime()


    '''
    ANALYSIS
    '''

    def analyze_experiment(self):
        # create lists for a ch1 sweep
        ch1_turns_sweep_list = []
        phonons_ch1_sweep_list = []

        # determine if a ch1 sweep occurred
        ch1_sweep_bool = len(self.phase_qvsa_ch1_turns_list) > 1
        rsb_sweep_bool = len(self.phase_qvsa_rsb_turns_list) > 1
        turns_sweep_bool = ch1_sweep_bool or rsb_sweep_bool

        # print results
        print("\tResults - QVSA Ramsey:")

        # sweep over ch1_turns
        for ch1_turns in self.phase_qvsa_ch1_turns_list:

            # handle errors from data processing
            try:
                # convert dataset to array
                dataset = np.reshape(self.results[np.where(self.results[:, 5] == ch1_turns), :],
                                     (-1, self.results.shape[1]))

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
                                                                                               self.repetitions,1)

                # get the phonons and instantiate the fitter class
                phonons = convert_ratios_to_coherent_phonons(ratios)
                fitter = fitSincGeneric()

                # format arguments for applet plotting
                if ch1_sweep_bool:
                    ccb_command = f'$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_heating_Ramsey_{ch1_turns}'
                    group = 'plotting.eggs_heating.ch1_sweep'
                    dataset_name = f'temp.plotting.results_eggs_heating_Ramsey_{ch1_turns}'
                    applet_name = f"EGGS Heating - Ramsey - CH1 Turns: {ch1_turns}"
                else:
                    ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_eggs_heating_Ramsey'
                    group = 'plotting.eggs_heating'
                    dataset_name = f'temp.plotting.results_eggs_heating_Ramsey'
                    applet_name = f"EGGS Heating - Ramsey"

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
                                    'rid': self.scheduler.rid}

                self.set_dataset('temp.plotting.results_eggs_heating_Ramsey_ch1_sweep',
                                 pyon.encode(plotting_results), broadcast=True)

                self.ccb.issue("create_applet", f"EGGS Heating - Ramsey",
                               '$python -m LAX_exp.applets.plot_matplotlib '
                               'temp.plotting.results_eggs_heating_Ramsey_ch1_sweep'
                               ' --num-subplots 1', group=['plotting', 'eggs_heating', 'ch1_sweep'])


            except Exception as e:
                print("Warning: unable to process data.")
                print(repr(e))