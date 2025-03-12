import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, SidebandCoolContinuous, QVSAPulse,
                                         QubitPulseShape, QubitRAP, Readout, RabiflopReadout, RescueIon)

from itertools import product
from artiq.coredevice import ad9910


class PuttermanPuzzle(LAXExperiment, Experiment):
    """
    Experiment: Putterman Puzzle

    Is a coherent state an eigenvalue of the annihilation operator?
    Let's find out.
    """
    name = 'Putterman Puzzle'
    kernel_invariants = {
        # subsequences etc.
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'rap_subsequence', 'rabiflop_subsequence',
        'profile_target_rap', 'profile_target_readout',

        # hardware values - core
        'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu', 'att_rap_mu',

        # hardware values - motional state prep

        # experiment/config related
        'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # get subsequences
        self.initialize_subsequence = InitializeQubit(self)
        self.sidebandcool_subsequence = SidebandCoolContinuous(self)
        self.readout_subsequence = Readout(self)

        # set target profile for stuff
        self.profile_target_rap =       5
        self.profile_target_readout =   6


        '''MOTIONAL STATE PREP - CONFIGURATION'''
        # QVSA pulse preparation
        self.motional_subsequence = QVSAPulse(self)

        '''RAP - CONFIGURATION'''
        self.setattr_argument("enable_rap",             BooleanValue(default=True), group='RAP')
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=101.3996, precision=6, step=1, min=50., max=400.), group="RAP")
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=10., precision=3, step=5, min=0.01, max=50), group="RAP")
        self.setattr_argument("time_rap_us",            NumberValue(default=90.5, precision=2, step=5, min=0.1, max=10000), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50), group="RAP")
        self.setattr_argument("att_rap_db",             NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5), group="RAP")

        '''HERALD - CONFIGURATION'''
        self.setattr_argument("enable_force_herald",    BooleanValue(default=True), group='herald')
        self.setattr_argument("force_herald_threshold", NumberValue(default=46, precision=0, step=10, min=0, max=10000), group='herald')

        '''RABIFLOP READOUT - CONFIGURATION'''
        self.rabiflop_subsequence = RabiflopReadout(self, profile_dds=self.profile_target_readout)

        # initialize other subsequences (which rely on build arguments-ish)
        self.rap_subsequence =  QubitRAP(self, ram_profile=self.profile_target_rap,
                                         ampl_max_pct=self.ampl_rap_pct, num_samples=500,
                                         pulse_shape="blackman")
        self.rescue_subsequence = RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CHECK INPUT FOR ERRORS
        '''
        # todo

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # RAP values
        self.freq_rap_center_ftw =  self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw =     self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu =          self.core.seconds_to_mu(self.time_rap_us * us)
        self.att_rap_mu =           att_to_mu(self.att_rap_db * dB)

        # heralding values
        self.time_force_herald_slack_mu = self.core.seconds_to_mu(150 * us)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            vals
            for vals in product(self.rabiflop_subsequence.time_readout_mu_list)
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)


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
        self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
        with self.core_dma.record('RAP_SUBSEQUENCE'):
            self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        # instantiate relevant variables
        counts_her = -1 # store heralded counts

        # retrieve relevant DMA sequences
        self.dma_handle_rap = self.core_dma.get_handle('RAP_SUBSEQUENCE')
        self.core.break_realtime()


        # MAIN EXECUTION LOOP
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                time_rabiflop_readout_mu =  config_vals[0]
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
                        self.qubit.set_att_mu(self.att_rap_mu)
                        self.core_dma.playback_handle(self.dma_handle_rap)

                    # herald ion via state-dependent fluorescence
                    self.readout_subsequence.run_dma()
                    self.pump.off()

                    # optional: force heralding (i.e. run until we succeed w/RAP)
                    if self.enable_force_herald:
                        counts_her = self.readout_subsequence.fetch_count()
                        if counts_her > self.force_herald_threshold:
                            self.core.break_realtime()
                            continue

                        # add minor slack if we proceed
                        at_mu(self.core.get_rtio_counter_mu() + self.time_force_herald_slack_mu)

                '''READOUT & STORE RESULTS'''
                # rabi flop & readout for motional detection
                self.rabiflop_subsequence.run_time(time_rabiflop_readout_mu)
                self.readout_subsequence.run_dma()

                # retrieve heralded measurement
                if not self.enable_force_herald:
                    counts_her = self.readout_subsequence.fetch_count()

                # retrieve readout results & update dataset
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(time_rabiflop_readout_mu,
                                    counts_res,
                                    counts_her)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        """
        Clean up the experiment.
        """
        self.core.break_realtime()


    '''
    HELPER FUNCTIONS
    '''
    # todo: not necessary now (but will be necessary if we wigner readout)


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

