import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, PHASE_MODE_TRACKING, _AD9910_REG_CFR1

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, RescueIon, SidebandCoolContinuous

from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper


class PSRSB(LAXExperiment, Experiment):
    """
    Experiment: Phase-Sensitive Red Sideband (PSRSB)

    Use QVSA to generate a phase-sensitive displacement, then apply
    the PSRSB technique for readout.
    """
    name = 'PSRSB'
    kernel_invariants = {
        # config/sweep
        'config_experiment_list', 'freq_psrsb_rsb_ftw_list', 'freq_psrsb_carrier_ftw_list',
        'phas_psrsb_carrier_pow_list',
        # QVSA/phaser related
        'freq_qvsa_carrier_hz', 'freq_qvsa_secular_hz', 'att_qvsa_mu', 'waveform_qvsa_pulseshape_vals',
        # PSRSB
        'profile_psrsb', 'att_qubit_mu', 'ampl_psrsb_rsb_asf', 'ampl_psrsb_carrier_asf', 'time_psrsb_rsb_mu',
        'time_psrsb_carrier_mu',
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'spinecho_wizard', 'pulse_shaper'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=100, precision=0, step=1, min=1, max=100000))

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandcool_subsequence =     SidebandCoolContinuous(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # QVSA configuration - pulse
        self.setattr_argument("freq_qvsa_carrier_mhz",      NumberValue(default=80., precision=6, step=10., min=0., max=1000.), group='QVSA')
        self.setattr_argument("freq_qvsa_secular_khz",      NumberValue(default=1281., precision=3, step=1., min=0., max=5000.), group='QVSA')
        self.setattr_argument("time_qvsa_us",               NumberValue(default=1000, precision=2, step=500, min=0.04, max=100000000), group='QVSA')
        self.setattr_argument("ampl_qvsa_pct_config",       PYONValue([1., 1., 1.]), group='QVSA', tooltip="[rsb_qvsa_pct, bsb_qvsa_pct, carrier_qvsa_pct]")
        self.setattr_argument("phase_qvsa_turns_config",    PYONValue([0., 0., 0.]), group='QVSA', tooltip="[rsb_qvsa_turns, bsb_qvsa_turns, carrier_qvsa_turns]")
        self.setattr_argument("att_qvsa_db",                NumberValue(default=31.5, precision=1, step=0.5, min=0, max=31.5), group='QVSA')

        # QVSA configuration - pulse shaping
        self.setattr_argument("enable_pulse_shaping",           BooleanValue(default=False), group='QVSA.pulse_shaping')
        self.setattr_argument("type_pulse_shape",               EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'), group='QVSA.pulse_shaping')
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000), group='QVSA.pulse_shaping')
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1000, precision=0, step=100, min=100, max=5000), group='QVSA.pulse_shaping')

        # PSRSB - RSB pulse
        self.setattr_argument("att_qubit_db",               NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5), group=self.name)
        self.setattr_argument("ampl_psrsb_rsb_pct",         NumberValue(default=50, precision=3, step=5, min=0.01, max=50), group=self.name)
        self.setattr_argument("time_psrsb_rsb_us",          NumberValue(default=148.57, precision=3, step=5, min=1, max=10000000), group=self.name)
        self.setattr_argument("freq_psrsb_rsb_mhz_list",    Scannable(
                                                                default=[
                                                                    CenterScan(100.7947, 0.01, 0.0001, randomize=True),
                                                                    ExplicitScan([100.7947]),
                                                                ],
                                                                global_min=60., global_max=200., global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group=self.name)

        self.setattr_argument("ampl_psrsb_carrier_pct",         NumberValue(default=50, precision=3, step=5, min=0.01, max=50), group=self.name)
        self.setattr_argument("time_psrsb_carrier_us",          NumberValue(default=8.51, precision=3, step=1, min=1, max=10000000), group=self.name)
        self.setattr_argument("freq_psrsb_carrier_mhz_list",    Scannable(
                                                                    default=[
                                                                        ExplicitScan([101.4181]),
                                                                        CenterScan(101.4181, 0.01, 0.0001, randomize=True),
                                                                    ],
                                                                    global_min=60., global_max=200., global_step=1,
                                                                    unit="MHz", scale=1, precision=6
                                                                ), group=self.name)
        self.setattr_argument("phas_psrsb_carrier_turns_list",  Scannable(
                                                                    default=[
                                                                        RangeScan(0, 1.0, 21, randomize=True),
                                                                        ExplicitScan([0.372]),
                                                                    ],
                                                                    global_min=-1.0, global_max=1.0, global_step=0.1,
                                                                    unit="turns", scale=1, precision=3
                                                                ), group=self.name)

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizard(self)
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        self.pulse_shaper = PhaserPulseShaper(self, np.array([0., 0., 0.5, 0., 0.]))

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE/VALIDATE INPUTS & CHECK ERRORS'''
        # ensure phaser oscillator amplitudes are configured correctly
        if (type(self.ampl_qvsa_pct_config) is not list) or (len(self.ampl_qvsa_pct_config) != 3):
            raise ValueError("Invalid QVSA oscillator amplitude configuration."
                             "Must be of list [rsb_ampl_pct, bsb_ampl_pct, carrier_ampl_pct].")
        elif not all(0. <= val <= 100. for val in self.ampl_qvsa_pct_config):
            raise ValueError("Invalid QVSA oscillator amplitude. Must be in range [0., 100.].")
        elif sum(self.ampl_qvsa_pct_config) >= 100.:
            raise ValueError("Invalid QVSA oscillator amplitudes. Total must sum to <= 100.")

        # ensure phaser oscillator phases are configured correctly
        if (type(self.phase_qvsa_turns_config) is not list) or (len(self.phase_qvsa_turns_config) != 3):
            raise ValueError("Invalid QVSA oscillator phase configuration."
                             "Must be of list [rsb_phas_turns, bsb_phas_turns, carrier_phas_turns].")

        # ensure phaser output frequency falls within valid DUC bandwidth
        if abs(self.freq_qvsa_carrier_mhz * MHz - self.phaser_eggs.freq_center_hz) >= 200. * MHz:
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        '''CONVERT VALUES TO MACHINE UNITS - QVSA'''
        # convert frequencies to Hz
        self.freq_qvsa_carrier_hz = self.freq_qvsa_carrier_mhz * MHz
        self.freq_qvsa_secular_hz = self.freq_qvsa_secular_khz * kHz

        # convert attenuation from dB to machine units
        self.att_qvsa_mu = att_to_mu(self.att_qvsa_db * dB)

        '''CONVERT VALUES TO MACHINE UNITS - PSRSB PULSES'''
        # convert attenuation from dB to machine units
        self.att_qubit_mu = att_to_mu(self.att_qubit_db * dB)

        # convert qubit DDS waveform values to machine units
        self.profile_psrsb = 6
        self.ampl_psrsb_rsb_asf =       self.qubit.amplitude_to_asf(self.ampl_psrsb_rsb_pct / 100.)
        self.ampl_psrsb_carrier_asf =   self.qubit.amplitude_to_asf(self.ampl_psrsb_carrier_pct / 100.)

        self.freq_psrsb_rsb_ftw_list =      np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                      for freq_mhz in list(self.freq_psrsb_rsb_mhz_list)])
        self.freq_psrsb_carrier_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                      for freq_mhz in list(self.freq_psrsb_carrier_mhz_list)])

        self.phas_psrsb_carrier_pow_list =  np.array([self.qubit.turns_to_pow(phas_turns)
                                                      for phas_turns in list(self.phas_psrsb_carrier_turns_list)])

        # convert time to machine units
        self.time_psrsb_rsb_mu =        self.core.seconds_to_mu(self.time_psrsb_rsb_us * us)
        self.time_psrsb_carrier_mu =    self.core.seconds_to_mu(self.time_psrsb_carrier_us * us)

        # create experiment config data structure
        self.config_experiment_list = np.stack(np.meshgrid(self.freq_psrsb_rsb_ftw_list,
                                                           self.freq_psrsb_carrier_ftw_list,
                                                           self.phas_psrsb_carrier_pow_list),
                                               -1, dtype=np.int32).reshape(-1, 3)
        np.random.shuffle(self.config_experiment_list)

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.waveform_qvsa_pulseshape_id = 0

        # set up the spin echo wizard generally
        # note: time_pulse_us is amount of time for each block
        self.spinecho_wizard.time_pulse_us =                self.time_qvsa_us
        self.spinecho_wizard.enable_pulse_shaping =         self.enable_pulse_shaping
        self.spinecho_wizard.pulse_shape_blocks =           False
        self.spinecho_wizard.type_pulse_shape =             self.type_pulse_shape
        self.spinecho_wizard.time_pulse_shape_rolloff_us =  self.time_pulse_shape_rolloff_us
        self.spinecho_wizard.freq_pulse_shape_sample_khz =  self.freq_pulse_shape_sample_khz
        self.spinecho_wizard.enable_delay_spinecho =        False
        self.spinecho_wizard.time_delay_spinecho_us =       250

        '''DESIGN WAVEFORM SEQUENCE'''
        # create bare waveform block sequence
        _sequence_blocks = np.zeros((1, 3, 2), dtype=float)

        # set oscillator amplitudes
        _sequence_blocks[:, 0, 0] = self.ampl_qvsa_pct_config[0]
        _sequence_blocks[:, 1, 0] = self.ampl_qvsa_pct_config[1]
        _sequence_blocks[:, 2, 0] = self.ampl_qvsa_pct_config[2]

        # set oscillator phases (accounting for oscillator delay time)
        phase_bsb_update_delay_turns = (self.freq_qvsa_secular_khz * kHz) * (self.phaser_eggs.t_sample_mu * ns)
        _sequence_blocks[:, 0, 1] = self.phase_qvsa_turns_config[0]
        _sequence_blocks[:, 1, 1] = self.phase_qvsa_turns_config[1] + phase_bsb_update_delay_turns
        _sequence_blocks[:, 2, 1] = self.phase_qvsa_turns_config[2]

        # create waveform
        self.spinecho_wizard.sequence_blocks = _sequence_blocks
        self.spinecho_wizard.calculate_pulseshape()
        self.spinecho_wizard.compile_waveform()

        # get waveform data and store in holding structure
        self.waveform_qvsa_pulseshape_vals = self.spinecho_wizard.get_waveform()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)


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
        # record phaser oscillator waveform
        # note: normally this would be encapsulated in phaser_record
        # e.g. in EGGSHeatingRDX-type exps
        self.core.break_realtime()
        delay_mu(1000000)  # add slack for recording DMA sequences (1000 us)
        self.waveform_qvsa_pulseshape_id = self.pulse_shaper.waveform_record(
            self.waveform_qvsa_pulseshape_vals[0],
            self.waveform_qvsa_pulseshape_vals[1],
            self.waveform_qvsa_pulseshape_vals[2]
        )
        self.core.break_realtime()

        # configure phaser carrier frequencies
        # self.phaser_configure(self.freq_qvsa_carrier_hz, self.freq_qvsa_secular_hz, self.phaser_eggs.phase_inherent_ch1_turns)
        self.phaser_eggs.frequency_configure(self.freq_qvsa_carrier_hz,
                                             [-self.freq_qvsa_secular_hz, self.freq_qvsa_secular_hz, 0., 0., 0.],
                                             self.phaser_eggs.phase_inherent_ch1_turns)
        self.core.break_realtime()

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
                freq_rsb_ftw =      config_vals[0]
                freq_carrier_ftw =  config_vals[1]
                phas_carrier_pow =  config_vals[2]
                self.core.break_realtime()

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                '''PHASE-SENSITIVE RED SIDEBAND SEQUENCE'''
                # create qvsa displacement
                self.phaser_run(self.waveform_qvsa_pulseshape_id)
                # run PSRSB detection
                self.psrsb_run(freq_rsb_ftw, freq_carrier_ftw, phas_carrier_pow)

                '''READOUT'''
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
                self.update_results(
                    freq_rsb_ftw,
                    counts,
                    freq_carrier_ftw,
                    phas_carrier_pow
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
    def psrsb_run(self, freq_rsb_ftw: TInt32, freq_carrier_ftw: TInt32, phas_carrier_pow: TInt32) -> TNone:
        """
        Run the phase-sensitive red-sideband detection sequence.
        Arguments:
            freq_rsb_ftw: RSB frequency in FTW.
            freq_carrier_ftw: Carrier frequency in FTW.
            phas_carrier_pow: Carrier phase (relative) in POW.
        """
        # set target profile and attenuation
        self.qubit.set_profile(self.profile_psrsb)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_qubit_mu)

        # synchronize start time to coarse RTIO clock
        time_start = now_mu() & ~0x7

        # run RSB pulse
        self.qubit.set_mu(freq_rsb_ftw, pow_=0, asf=self.ampl_psrsb_rsb_asf,
                          profile=self.profile_psrsb,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_start)
        self.qubit.on()
        delay_mu(self.time_psrsb_rsb_mu)
        self.qubit.off()

        # run carrier pulse
        self.qubit.set_mu(freq_carrier_ftw, pow_=phas_carrier_pow, asf=self.ampl_psrsb_carrier_asf,
                          profile=self.profile_psrsb,
                          phase_mode=PHASE_MODE_TRACKING, ref_time_mu=time_start)
        self.qubit.on()
        delay_mu(self.time_psrsb_carrier_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main EGGS pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_qvsa_mu, self.att_qvsa_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # EGGS - STOP
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        todo: document
        """
        # todo: implement
        pass
