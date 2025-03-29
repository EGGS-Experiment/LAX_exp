import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon,
    SidebandCoolContinuous, SidebandReadout
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper

# todo: bring back ampl sweep


class DynamicAntisqueezing(LAXExperiment, Experiment):
    """
    Experiment: Dynamic Antisqueezing

    Antisqueeze using a different waveform than the initial squeezing.
    """
    name = 'Dynamic Antisqueezing'
    kernel_invariants = {
        # config/sweep
        'config_experiment_list', 'freq_sideband_readout_ftw_list', 'time_readout_mu_list',
        'freq_pulse_secular_hz_list',

        # squeezing
        'freq_squeeze_carrier_hz', 'att_squeeze_mu', 'waveform_squeezing_compiled',

        # antisqueezing
        'freq_antisqueeze_carrier_hz', 'att_antisqueeze_mu', 'phase_antisqueeze_offset_turns_list',
        'waveform_idx_to_phase_turns', 'waveform_idx_to_compiled',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'sidebandreadout_subsequence', 'rescue_subsequence', 'spinecho_wizard', 'pulse_shaper'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=10000, precision=0, step=1, min=1, max=100000))

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandcool_subsequence =     SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # readout
        self.setattr_argument("time_readout_us_list",   Scannable(
                                                            default=[
                                                                ExplicitScan([120.5]),
                                                                RangeScan(0, 1500, 100, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group='sideband_readout')

        # general pulse configuration
        self.setattr_argument("freq_pulse_secular_khz_list",    Scannable(
                                                                    default=[
                                                                        CenterScan(777.5, 4, 0.5, randomize=True),
                                                                        ExplicitScan([1276.15]),
                                                                        ExplicitScan([767.2, 319.2, 1582, 3182]),
                                                                    ],
                                                                    global_min=0., global_max=10000, global_step=1,
                                                                    unit="kHz", scale=1, precision=3
                                                                ), group='pulse.general')
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=True), group='pulse.general')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'), group='pulse.general')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000), group='pulse.general')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1000, precision=0, step=100, min=100, max=5000), group='pulse.general')

        # pulse 0 (squeezing) - configuration
        self.setattr_argument("enable_squeezing",           BooleanValue(default=True), group='pulse.squeeze')
        self.setattr_argument("freq_squeeze_carrier_mhz",   NumberValue(default=80., precision=6, step=10., min=0., max=1000.), group='pulse.squeeze')
        self.setattr_argument("att_squeeze_db",             NumberValue(default=10., precision=1, step=0.5, min=0, max=31.5), group='pulse.squeeze')
        self.setattr_argument("time_squeeze_us",            NumberValue(default=200, precision=2, step=500, min=0.04, max=100000000), group='pulse.squeeze')
        self.setattr_argument("ampl_squeeze_pct_config",    PYONValue([40., 40.]), group='pulse.squeeze', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("phase_squeeze_turns_config", PYONValue([0.0, 0.0]), group='pulse.squeeze', tooltip="[rsb_turns, bsb_turns]")

        # pulse 1 (antisqueezing) - configuration
        self.setattr_argument("enable_antisqueezing",       BooleanValue(default=True), group='pulse.antisqueeze')
        self.setattr_argument("freq_antisqueeze_carrier_mhz",   NumberValue(default=40., precision=6, step=10., min=0., max=1000.), group='pulse.antisqueeze')
        self.setattr_argument("att_antisqueeze_db",         NumberValue(default=4., precision=1, step=0.5, min=0, max=31.5), group='pulse.antisqueeze')
        self.setattr_argument("time_antisqueeze_us",        NumberValue(default=100, precision=2, step=500, min=0.04, max=100000000), group='pulse.antisqueeze')
        self.setattr_argument("ampl_antisqueeze_pct_config",    PYONValue([45., 45.]), group='pulse.antisqueeze', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("phase_antisqueeze_turns_config", PYONValue([0.0, 0.5]), group='pulse.antisqueeze', tooltip="[rsb_turns, bsb_turns]")

        # pulse 1 (antisqueezing) - sweep
        # self.setattr_argument("ampl_antisqueeze_scaling_frac_list", Scannable(
        #                                                                 default=[
        #                                                                     ExplicitScan([1.0]),
        #                                                                     RangeScan(0., 1.0, 21, randomize=True),
        #                                                                 ],
        #                                                                 global_min=0., global_max=10.0, global_step=0.1,
        #                                                                 unit="frac", scale=1, precision=4
        #                                                             ), group='pulse.antisqueeze.sweep')
        self.setattr_argument("target_antisqueeze_phase",               EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB'), group='pulse.antisqueeze.sweep')
        self.setattr_argument("phase_antisqueeze_offset_turns_list",    Scannable(
                                                                            default=[
                                                                                RangeScan(0, 1.0, 21, randomize=True),
                                                                                ExplicitScan([0.372]),
                                                                            ],
                                                                            global_min=-1.0, global_max=1.0, global_step=0.1,
                                                                            unit="turns", scale=1, precision=3
                                                                        ), group='pulse.antisqueeze.sweep')

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)
        # set correct phase delays for field geometries
        self.pulse_shaper = PhaserPulseShaper(self, np.array([0., 0., 0., 0., 0.]))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        # check experiment arguments for errors
        self._prepare_argument_checks()

        '''CONVERT VALUES TO MACHINE UNITS - 729nm READOUT'''
        self.freq_sideband_readout_ftw_list =   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.time_readout_mu_list =             np.array([self.core.seconds_to_mu(time_us * us)
                                                          for time_us in self.time_readout_us_list])

        '''CONVERT VALUES TO MACHINE UNITS - SQUEEZING/ANTISQUEEZING'''
        # convert frequencies to Hz
        self.freq_squeeze_carrier_hz = self.freq_squeeze_carrier_mhz * MHz
        self.freq_antisqueeze_carrier_hz = self.freq_antisqueeze_carrier_mhz * MHz
        self.freq_pulse_secular_hz_list = np.array(list(self.freq_pulse_secular_khz_list)) * kHz

        # convert attenuation from dB to machine units
        self.att_squeeze_mu = att_to_mu(self.att_squeeze_db * dB)
        self.att_antisqueeze_mu = att_to_mu(self.att_antisqueeze_db * dB)

        '''FINISH EXPERIMENT CONFIG'''
        # create empty holder for compiled squeezing waveform
        self.waveform_squeezing_compiled = None
        self.waveform_squeezing_id = -1

        # map antisqueezing waveform to array index for programmatic DMA recording/playback
        self.phase_antisqueeze_offset_turns_list = np.array(list(self.phase_antisqueeze_offset_turns_list))
        # used to map an index to the phase value (-ish)
        self.waveform_idx_to_phase_turns =      np.arange(len(self.phase_antisqueeze_offset_turns_list))
        self.waveform_idx_to_compiled =         list()      # store compiled waveforms
        self.waveform_idx_to_pulseshaper_id =   np.zeros(len(self.phase_antisqueeze_offset_turns_list),
                                                         dtype=np.int32) # store pulse shaper waveform ID

        # create experiment config data structure
        self.config_experiment_list = np.stack(np.meshgrid(self.freq_sideband_readout_ftw_list,
                                                           self.freq_pulse_secular_hz_list,
                                                           self.waveform_idx_to_phase_turns,
                                                           self.time_readout_mu_list),
                                               -1, dtype=float).reshape(-1, 4)
        np.random.shuffle(self.config_experiment_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        '''SANITIZE/VALIDATE INPUTS & CHECK ERRORS'''
        if self.enable_squeezing:
            # ensure phaser output frequency falls within valid DUC bandwidth - squeezing
            if abs(self.freq_squeeze_carrier_mhz * MHz - self.phaser_eggs.freq_center_hz) >= 200. * MHz:
                raise ValueError("Error: squeezing frequency outside +/- 300 MHz phaser DUC bandwidth.")

            # ensure phaser oscillator amplitudes are configured correctly - squeezing
            if (type(self.ampl_squeeze_pct_config) is not list) or (len(self.ampl_squeeze_pct_config) != 2):
                raise ValueError("Invalid squeezing oscillator amplitude configuration."
                                 "Must be of list [rsb_ampl_pct, bsb_ampl_pct].")
            elif not all(0. <= val <= 100. for val in self.ampl_squeeze_pct_config):
                raise ValueError("Invalid squeezing oscillator amplitude. Must be in range [0., 100.].")
            elif sum(self.ampl_squeeze_pct_config) >= 100.:
                raise ValueError("Invalid squeezing oscillator amplitudes. Total must sum to <= 100.")

            # ensure phaser oscillator phases are configured correctly - squeezing
            if (type(self.phase_squeeze_turns_config) is not list) or (len(self.phase_squeeze_turns_config) != 2):
                raise ValueError("Invalid squeezing oscillator phase configuration."
                                 "Must be of list [rsb_phas_turns, bsb_phas_turns].")

        if self.enable_antisqueezing:
            # ensure phaser output frequency falls within valid DUC bandwidth - antisqueezing
            if abs(self.freq_antisqueeze_carrier_mhz * MHz - self.phaser_eggs.freq_center_hz) >= 200. * MHz:
                raise ValueError("Error: antisqueezing frequency outside +/- 300 MHz phaser DUC bandwidth.")

            # ensure phaser oscillator amplitudes are configured correctly - antisqueezing
            if (type(self.ampl_antisqueeze_pct_config) is not list) or (len(self.ampl_antisqueeze_pct_config) != 2):
                raise ValueError("Invalid antisqueezing oscillator amplitude configuration."
                                 "Must be of list [rsb_ampl_pct, bsb_ampl_pct].")
            elif not all(0. <= val <= 100. for val in self.ampl_antisqueeze_pct_config):
                raise ValueError("Invalid antisqueezing oscillator amplitude. Must be in range [0., 100.].")
            elif sum(self.ampl_antisqueeze_pct_config) >= 100.:
                raise ValueError("Invalid antisqueezing oscillator amplitudes. Total must sum to <= 100.")
            # tmp remove: for sweeping ampl scaling
            # elif sum(self.ampl_antisqueeze_pct_config) * max(1., max(list(self.ampl_antisqueeze_scaling_frac_list))) >= 100.:
            #     raise ValueError("Invalid antisqueezing oscillator amplitudes. Scaled total must sum to <= 100.")
            # tmp remove: for sweeping ampl scaling

            # ensure phaser oscillator phases are configured correctly - squeezing
            if (type(self.phase_antisqueeze_turns_config) is not list) or (len(self.phase_antisqueeze_turns_config) != 2):
                raise ValueError("Invalid antisqueezing oscillator phase configuration."
                                 "Must be of list [rsb_phas_turns, bsb_phas_turns].")

        # tmp remove: for sweeping ampl scaling
        # # ensure only 1 waveform sweep at a time to reduce DMA overhead
        # if (len(self.ampl_antisqueeze_scaling_frac_list) != 1) and (len(self.phase_antisqueeze_offset_turns_list) != 1):
        #     raise ValueError("Only a single waveform sweep allowed at a time to reduce DMA overhead.")
        # tmp remove: for sweeping ampl scaling

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE SQUEEZING WAVEFORM'''
        _sequence_block_squeeze = [
            {
                "oscillator_parameters": [
                    [self.ampl_squeeze_pct_config[0], self.phase_squeeze_turns_config[0]],
                    [self.ampl_squeeze_pct_config[1], self.phase_squeeze_turns_config[1]],
                ],
                "config": {
                    "time_us": self.time_squeeze_us,
                    "pulse_shaping": True,
                    "pulse_shaping_config": {
                        "pulse_shape":          "sine_squared",
                        "pulse_shape_rising":   self.enable_pulse_shaping,
                        "pulse_shape_falling":  self.enable_pulse_shaping,
                        "sample_rate_khz":      self.freq_pulse_shape_sample_khz,
                        "rolloff_time_us":      self.time_pulse_shape_rolloff_us
                    }
                }
            },
        ]

        # selectively disable by using 0 amplitude
        if not self.enable_squeezing:
            _sequence_block_squeeze[0]["oscillator_parameters"][0][0] = 0.
            _sequence_block_squeeze[0]["oscillator_parameters"][1][0] = 0.

        # create squeezing waveform
        self.spinecho_wizard.sequence_blocks = _sequence_block_squeeze
        self.spinecho_wizard.compile_waveform()

        # get waveform data and store in holding structure
        self.waveform_squeezing_compiled = self.spinecho_wizard.get_waveform()


        '''PREPARE WAVEFORM COMPILATION - ANTISQUEEZING'''
        _sequence_block_antisqueeze = [
            {
                "oscillator_parameters": [
                    [self.ampl_antisqueeze_pct_config[0], self.phase_antisqueeze_turns_config[0]],
                    [self.ampl_antisqueeze_pct_config[1], self.phase_antisqueeze_turns_config[1]],
                ],
                "config": {
                    "time_us": self.time_antisqueeze_us,
                    "pulse_shaping": self.type_pulse_shape,
                    "pulse_shaping_config": {
                        "pulse_shape": "sine_squared",
                        "pulse_shape_rising":   self.enable_pulse_shaping,
                        "pulse_shape_falling":  self.enable_pulse_shaping,
                        "sample_rate_khz":      self.freq_pulse_shape_sample_khz,
                        "rolloff_time_us":      self.time_pulse_shape_rolloff_us
                    }
                }
            },
        ]
        
        '''DESIGN WAVEFORM SEQUENCE'''
        # selectively disable by using 0 amplitude
        if not self.enable_antisqueezing:
            _sequence_block_antisqueeze[0]["oscillator_parameters"][0][0] = 0.
            _sequence_block_antisqueeze[0]["oscillator_parameters"][1][0] = 0.

        # adjust oscillator phases based on user configuration
        if self.target_antisqueeze_phase == "RSB":
            phas_update_arr = [1., 0.]
        elif self.target_antisqueeze_phase == "BSB":
            phas_update_arr = [0., 1.]
        elif self.target_antisqueeze_phase == "RSB-BSB":
            phas_update_arr = [1., -1.]
        elif self.target_antisqueeze_phase == "RSB+BSB":
            phas_update_arr = [1., 1.]

        # record EGGS pulse waveforms
        for i, phas_turns in enumerate(self.phase_antisqueeze_offset_turns_list):
            # update sequence block with target phase
            _sequence_block_antisqueeze[0]["oscillator_parameters"][0][1] += phas_update_arr[0] * phas_turns
            _sequence_block_antisqueeze[0]["oscillator_parameters"][1][1] += phas_update_arr[1] * phas_turns

            # create waveform
            self.spinecho_wizard.sequence_blocks = _sequence_block_antisqueeze
            self.spinecho_wizard.compile_waveform()

            # get waveform data and store in holding structure
            self.waveform_idx_to_compiled.append(self.spinecho_wizard.get_waveform())

        # tmp remove
        # _wav_print_idk = self.waveform_idx_to_compiled[-1]
        # print(_wav_print_idk[0])
        # print(_wav_print_idk[1])
        # print(_wav_print_idk[2])
        # print(_sequence_blocks)
        # self.spinecho_wizard.display_waveform()
        # self.set_dataset("waveforms", self.waveform_idx_to_compiled)
        # tmp remove

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    '''
    MAIN SEQUENCE
    '''
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
        self.core.break_realtime()

        # load waveform DMA handles
        self.pulse_shaper.waveform_load()
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                freq_sideband_hz =  config_vals[1]
                phas_wav_idx =      np.int32(config_vals[2])
                time_readout_mu =   np.int64(config_vals[3])
                self.core.break_realtime()

                # get corresponding RSB phase and waveform ID from the index
                phas_offset_turns = self.phase_antisqueeze_offset_turns_list[phas_wav_idx]
                waveform_id = self.waveform_idx_to_pulseshaper_id[phas_wav_idx]
                self.core.break_realtime()

                # configure EGGS tones and set readout frequency
                self.phaser_eggs.frequency_configure(self.freq_squeeze_carrier_hz,
                                                     [-freq_sideband_hz, freq_sideband_hz, 0., 0., 0.],
                                                     self.phaser_eggs.phase_inherent_ch1_turns)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                '''SQUEEZE + ANTISQUEEZE'''
                self.phaser_run(waveform_id)

                '''READOUT'''
                self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
                self.update_results(
                    freq_readout_ftw,
                    counts,
                    freq_sideband_hz,
                    phas_offset_turns,
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

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main EGGS pulse together with supporting functionality.
        Arguments:
            waveform_id: the ID of the waveform to run.
        """
        '''START/SETUP'''
        self.phaser_eggs.phaser_setup(self.att_squeeze_mu, self.att_squeeze_mu)
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()

        '''SQUEEZE'''
        self.pulse_shaper.waveform_playback(self.waveform_squeezing_id)

        '''ANTISQUEEZE'''
        # reconfigure attenuation for antisqueezing
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(self.att_antisqueeze_mu)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(self.att_antisqueeze_mu)

        # reconfigure carrier freq for antisqueezing
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(self.freq_antisqueeze_carrier_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(self.freq_antisqueeze_carrier_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.duc_stb()

        # reconfigure CH1 phase for antisqueezing
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(self.phaser_eggs.phase_inherent_ch1_turns +
                                                  (self.freq_antisqueeze_carrier_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns))
        self.phaser_eggs.duc_stb()

        # run antisqueezing pulse
        self.pulse_shaper.waveform_playback(waveform_id)

        '''STOP'''
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record squeezing waveform onto DMA
        _wav_data_ampl, _wav_data_phas, _wav_data_time = self.waveform_squeezing_compiled
        self.core.break_realtime()
        delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
        self.waveform_squeezing_id = self.pulse_shaper.waveform_record(_wav_data_ampl,
                                                                       _wav_data_phas,
                                                                       _wav_data_time)
        self.core.break_realtime()

        # record antisqueezing waveforms onto DMA
        for i in range(len(self.phase_antisqueeze_offset_turns_list)):
            # get waveform for given RSB phase
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self.waveform_idx_to_compiled[i]
            self.core.break_realtime()

            # record phaser pulse sequence and save returned waveform ID
            delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
            self.waveform_idx_to_pulseshaper_id[i] = self.pulse_shaper.waveform_record(_wav_data_ampl,
                                                                                       _wav_data_phas,
                                                                                       _wav_data_time)
            self.core.break_realtime()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        todo: document
        """
        # todo: implement
        pass
