import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, QVSAPulse, QubitRAP, Readout,
    ReadoutAdaptive, RabiflopReadout, RescueIon
)

from itertools import product


class PuttermanPuzzle(LAXExperiment, Experiment):
    """
    Experiment: Putterman Puzzle

    Is a coherent state an eigenvalue of the annihilation operator?
    Let's find out.
    """
    name = 'Putterman Puzzle'
    kernel_invariants = {
        # subsequences etc.
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'readout_adaptive_subsequence',
        'rescue_subsequence', 'rap_subsequence', 'rabiflop_subsequence',

        # hardware values - core
        'freq_rap_dev_ftw', 'time_rap_mu', 'att_rap_mu', 'time_force_herald_slack_mu',

        # experiment/config related
        'profile_729_SBC', 'profile_729_rap', 'profile_729_readout', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=10, precision=0, step=1, min=1, max=100000))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_SBC =      4
        self.profile_729_rap =      5
        self.profile_729_readout =  6

        '''MOTIONAL STATE PREP - CONFIGURATION'''
        # ram-based continuous sideband cooling
        self.sidebandcool_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )

        # QVSA pulse
        self.motional_subsequence = QVSAPulse(self)

        '''RAP - CONFIGURATION'''
        self.setattr_argument("enable_rap",             BooleanValue(default=True), group='RAP')
        self.setattr_argument("enable_quench",          BooleanValue(default=True), group='RAP')
        self.setattr_argument("freq_rap_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.3977]),
                                                                    CenterScan(101.3977, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.9801, 101.9901, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="RAP")
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=12.5, precision=3, step=5, min=0.01, max=1000), group="RAP")
        self.setattr_argument("time_rap_us",            NumberValue(default=125, precision=2, step=5, min=0.1, max=10000), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50), group="RAP")
        self.setattr_argument("att_rap_db",             NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5), group="RAP")

        '''HERALD - CONFIGURATION'''
        self.setattr_argument("enable_herald",          BooleanValue(default=True), group='herald')

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
        self.readout_subsequence =      Readout(self)
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
        freq_rap_center_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_rap_center_mhz_list], dtype=np.int32)
        self.freq_rap_dev_ftw =     self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu =          self.core.seconds_to_mu(self.time_rap_us * us)
        self.att_rap_mu =           att_to_mu(self.att_rap_db * dB)

        # heralding/readout
        self.time_force_herald_slack_mu = self.core.seconds_to_mu(150 * us)
        freq_rabiflop_readout_ftw_list =    np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                      for freq_mhz in self.freq_rabiflop_readout_mhz_list],
                                                     dtype=np.int32)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            vals
            for vals in product(
                freq_rap_center_ftw_list,
                self.rabiflop_subsequence.time_readout_mu_list,
                freq_rabiflop_readout_ftw_list
            )
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

        # # tmp remove - high fock test
        self.freq_bsb_ftw =     self.qubit.frequency_to_ftw(101.1065 * MHz)
        self.freq_rsb_ftw =     self.qubit.frequency_to_ftw(100.7554 * MHz)
        self.ampl_fock_asf =    self.qubit.amplitude_to_asf(0.5)
        self.att_fock_mu =      att_to_mu(8. * dB)
        self.time_fock_mu =     self.core.seconds_to_mu(2.5 * us)
        self.profile_fock =     3
        # # tmp remove - high fock test

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

        # configure RAP & record onto DMA
        with self.core_dma.record('RAP_SUBSEQUENCE'):
            # ensure correct att set (RAP subseq doesn't do this for us)
            self.qubit.set_att_mu(self.att_rap_mu)
            self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        # instantiate relevant variables
        ion_state = (-1, 0, np.int64(0))

        # retrieve relevant DMA sequences.handles
        self.motional_subsequence.pulse_shaper.waveform_load()
        self.core.break_realtime()
        dma_handle_rap = self.core_dma.get_handle('RAP_SUBSEQUENCE')
        self.core.break_realtime()


        # MAIN EXECUTION LOOP
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_rap_center_ftw =       np.int32(config_vals[0])
                time_rabiflop_readout_mu =  config_vals[1]
                freq_rabiflop_readout_ftw = np.int32(config_vals[2])
                self.core.break_realtime()

                # configure RAP pulse
                self.rap_subsequence.configure(self.time_rap_mu, freq_rap_center_ftw, self.freq_rap_dev_ftw)
                self.core.break_realtime()

                # configure rabiflop readout frequency
                self.qubit.set_mu(freq_rabiflop_readout_ftw, asf=self.qubit.ampl_qubit_asf,
                                  profile=self.profile_729_readout)
                self.core.break_realtime()

                '''PREPARE TARGET MOTIONAL STATE'''
                while True:

                    '''INITIALIZE'''
                    # initialize ion in S-1/2 state, then SBC to ground motional state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()

                    '''APPLY MOTIONAL INTERACTION'''
                    # use QVSA for motional excitation
                    self.motional_subsequence.run_pulse()

                    '''APPLY a^\dag (VIA RAP + HERALD)'''
                    if self.enable_rap:
                        # run RAP on motional sideband (a^\dag operator)
                        self.core_dma.playback_handle(dma_handle_rap)

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
                # configurable quench: put spin state back into ground state
                if self.enable_quench:
                    # note: ensure DDS set to readout parameters; should already be
                    # in correct profile b/c we do heralding immediately before
                    self.pump.readout()
                    self.repump_qubit.on()
                    delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                    self.repump_qubit.off()

                # rabi flop & readout for motional detection
                self.rabiflop_subsequence.run_time(time_rabiflop_readout_mu)
                self.readout_subsequence.run_dma()

                # retrieve readout results & update dataset
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(freq_rap_center_ftw,
                                    time_rabiflop_readout_mu,
                                    freq_rabiflop_readout_ftw,
                                    counts_res)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_fock(self, freq_ftw: TInt32, time_mu: TInt64) -> TNone:
        # with parallel:
        #     self.pump.readout()
        #     with sequential:
        #         self.qubit.set_mu(freq_ftw, asf=self.ampl_fock_asf, profile=self.profile_fock)
        #         self.qubit.set_profile(self.profile_fock)
        #         self.qubit.io_update()
        #         self.qubit.set_att_mu(self.att_fock_mu)

        self.qubit.set_mu(freq_ftw, asf=self.ampl_fock_asf, profile=self.profile_fock)
        self.qubit.set_profile(self.profile_fock)
        self.qubit.io_update()
        self.qubit.set_att_mu(self.att_fock_mu)

        self.qubit.on()
        delay_mu(time_mu)
        self.qubit.off()

        # self.repump_qubit.on()
        # delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
        # self.repump_qubit.off()

