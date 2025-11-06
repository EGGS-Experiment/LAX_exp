from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from numpy import int32, int64
from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, QVSAPulse, QubitRAP,
    ReadoutAdaptive, RabiflopReadout, RescueIon
)

# todo: make fast RSB a real configurable option


class PuttermanPuzzle(LAXExperiment, Experiment):
    """
    Experiment: Putterman Puzzle

    Is a coherent state an eigenvalue of the annihilation operator?
    Let's find out.
    """
    name = 'Putterman Puzzle'
    kernel_invariants = {
        # subsequences etc.
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_adaptive_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'rabiflop_subsequence',

        # hardware values - core
        'freq_rap_dev_ftw', 'time_rap_mu', 'att_rap_mu', 'time_force_herald_slack_mu',

        # experiment/config related
        'profile_729_frsb', 'profile_729_SBC', 'profile_729_rap', 'profile_729_readout', 'config_experiment_list',

        # tmp remove - frsb
        'time_frsb_mu', 'ampl_rap_asf'
        # tmp remove - frsb
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=75, precision=0, step=1, min=1, max=100000))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_frsb =     3
        self.profile_729_SBC =      4
        self.profile_729_rap =      5
        self.profile_729_readout =  6

        '''MOTIONAL STATE PREP - CONFIGURATION'''
        # ram-based continuous sideband cooling
        self.sidebandcool_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        # QVSA pulse
        self.motional_subsequence = QVSAPulse(self)

        '''A OPERATOR - CONFIGURATION'''
        self.setattr_argument("enable_rap",     BooleanValue(default=True), group='a_hat')
        self.setattr_argument("enable_quench",  BooleanValue(default=True), group='a_hat')
        self.setattr_argument("enable_herald",  BooleanValue(default=True), group='a_hat')
        self.setattr_argument("freq_rap_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.1012]),
                                                                    CenterScan(101.1012, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.0912, 101.1112, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="a_hat")
        self.setattr_argument("freq_rap_dev_khz",   NumberValue(default=175, precision=3, step=5, min=0.01, max=1000, unit="kHz", scale=1.),
                              group="a_hat")
        self.setattr_argument("time_rap_us",    NumberValue(default=1000, precision=2, step=5, min=0.1, max=100000, unit="us", scale=1.),
                              group="a_hat")
        self.setattr_argument("ampl_rap_pct",   NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="a_hat")
        self.setattr_argument("att_rap_db",     NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.),
                              group="a_hat")

        '''RABIFLOP READOUT - CONFIGURATION'''
        self.setattr_argument("freq_rabiflop_readout_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.9851]),
                                                                    CenterScan(101.9851, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.9801, 101.9901, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="rabiflop_readout")
        self.rabiflop_subsequence = RabiflopReadout(self, profile_dds=self.profile_729_readout)

        # initialize all other subsequences (which rely on build
        # arguments-ish, or don't themselves have arguments)
        self.rap_subsequence =  QubitRAP(
            self, ram_profile=self.profile_729_rap, ram_addr_start=202, num_samples=500,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_adaptive_subsequence = ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence =       RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # RAP values
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.ampl_rap_asf =     self.qubit.amplitude_to_asf(self.ampl_rap_pct / 100.)
        self.time_rap_mu =      self.core.seconds_to_mu(self.time_rap_us * us)
        self.att_rap_mu =       self.qubit.cpld.att_to_mu(self.att_rap_db * dB)
        freq_rap_center_ftw_list =  [self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                     for freq_mhz in self.freq_rap_center_mhz_list]

        # heralding/readout
        self.time_force_herald_slack_mu =   self.core.seconds_to_mu(150 * us)
        freq_rabiflop_readout_ftw_list =    [self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                             for freq_mhz in self.freq_rabiflop_readout_mhz_list]


        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create experiment config
        self.config_experiment_list = create_experiment_config(
            freq_rap_center_ftw_list,
            self.rabiflop_subsequence.time_readout_mu_list,
            freq_rabiflop_readout_ftw_list,
            config_type=int64,
            shuffle_config=True
        )

        # tmp remove - fast RSB
        time_frsb_us = 6.22
        self.set_dataset("time_frsb_us", time_frsb_us)
        self.time_frsb_mu = self.core.seconds_to_mu(time_frsb_us * us)
        # tmp remove - fast RSB

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

        # configure RAP & record onto DMA
        with self.core_dma.record('RAP_SUBSEQUENCE'):
            # ensure correct att set (RAP subseq doesn't do this for us)
            self.qubit.set_att_mu(self.att_rap_mu)
            self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        ion_state = (-1, 0, int64(0)) # store ion state from ReadoutAdaptive

        # retrieve relevant DMA sequences.handles
        self.motional_subsequence.pulse_shaper.waveform_load()
        self.core.break_realtime()
        dma_handle_rap = self.core_dma.get_handle('RAP_SUBSEQUENCE')

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_rap_center_ftw =       int32(config_vals[0])
                time_rabiflop_readout_mu =  config_vals[1]
                freq_rabiflop_readout_ftw = int32(config_vals[2])
                self.core.break_realtime()

                # configure adag pulse
                if self.enable_rap:
                    self.rap_subsequence.configure(self.time_rap_mu, freq_rap_center_ftw, self.freq_rap_dev_ftw)
                else:
                    self.qubit.set_mu(freq_rap_center_ftw, asf=self.ampl_rap_asf,
                                      profile=self.profile_729_frsb,
                                      phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(50000)
                # configure rabiflop readout frequency
                self.qubit.set_mu(freq_rabiflop_readout_ftw, asf=self.qubit.ampl_qubit_asf,
                                  profile=self.profile_729_readout,
                                  phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(25000)

                '''PREPARE TARGET MOTIONAL STATE'''
                while True:

                    '''INITIALIZE MOTIONAL STATE'''
                    # initialize ion in S-1/2 state, then SBC to ground motional state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()
                    # use QVSA to generate motional excitation
                    self.motional_subsequence.run_pulse()

                    '''APPLY a^\dag (VIA RAP + HERALD)'''
                    if self.enable_rap:
                        self.core_dma.playback_handle(dma_handle_rap)
                    else:
                        self.pulse_fock(self.time_frsb_mu)

                    # optional: herald ion via state-dependent fluorescence
                    if self.enable_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(20000)
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0: continue
                        # otherwise, add minor slack and proceed
                        at_mu(self.core.get_rtio_counter_mu() + self.time_force_herald_slack_mu)

                    # force break loop by default
                    break

                '''READOUT & STORE RESULTS'''
                # configurable quench: quench spin back to S-1/2 for characteristic readout
                if self.enable_quench:
                    self.initialize_subsequence.quench()

                # rabi flop & readout for motional detection
                self.rabiflop_subsequence.run_time(time_rabiflop_readout_mu)

                # retrieve readout results & update dataset
                ion_state = self.readout_adaptive_subsequence.run()
                self.update_results(freq_rap_center_ftw,
                                    ion_state[0],
                                    time_rabiflop_readout_mu,
                                    freq_rabiflop_readout_ftw)
                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_fock(self, time_mu: TInt64) -> TNone:
        self.qubit.set_profile(self.profile_729_frsb)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_rap_mu)
        delay_mu(5000)
        self.qubit.on()
        delay_mu(time_mu)
        self.qubit.off()

