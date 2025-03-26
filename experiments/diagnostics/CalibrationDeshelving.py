import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, QubitRAP, Readout, RescueIon
# todo: ampl calib???


class CalibrationDeshelving(LAXExperiment, Experiment):
    """
    Calibration: Deshelving

    Calibrate/measure deshelving (i.e. 854nm) pulse parameters.
    """
    name = 'Calibration Deshelving'
    kernel_invariants = {
        # hardware objects
        'initialize_subsequence', 'rap_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # hardware parameters
        'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu', 'att_rap_mu',
        'time_force_herald_slack_mu',

        # experiment values
        'profile_854_target', 'profile_729_rap', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=40, precision=0, step=1, min=1, max=100000))

        '''DESHELVING - CONFIGURATION'''
        # allocate DDS profiles
        self.profile_854_target =   6
        self.profile_729_rap =      1

        # deshelving scan parameters
        self.setattr_argument("freq_deshelve_mhz_list", Scannable(
                                                        default=[
                                                            CenterScan(110, 5, 0.2, randomize=True),
                                                            ExplicitScan([110.]),
                                                            RangeScan(90., 130., 80, randomize=True),
                                                        ],
                                                        global_min=60, global_max=200, global_step=0.01,
                                                        unit="MHz", scale=1, precision=6
                                                    ), group=self.name)
        self.setattr_argument("ampl_deshelve_pct_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([20]),
                                                            RangeScan(1, 50, 50, randomize=True),
                                                        ],
                                                        global_min=0.1, global_max=50, global_step=10,
                                                        unit="pct", scale=1, precision=3
                                                    ), group=self.name)
        self.setattr_argument("time_deshelve_us_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([20]),
                                                            RangeScan(10, 100, 91, randomize=True),
                                                        ],
                                                        global_min=1, global_max=2000, global_step=10,
                                                        unit="us", scale=1, precision=3
                                                    ), group=self.name)

        '''STATE PREPARATION - CONFIGURATION'''
        # RAP parameters: used to prepare D-5/2 dark state
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=101.3977, precision=6, step=1, min=50., max=400.), group="RAP")
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=12.5, precision=3, step=5, min=0.01, max=1000), group="RAP")
        self.setattr_argument("time_rap_us",            NumberValue(default=125, precision=2, step=5, min=0.1, max=10000), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50), group="RAP")
        self.setattr_argument("att_rap_db",             NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5), group="RAP")

        # heralding parameters: used to ensure successful state preparation
        self.setattr_argument("enable_herald",          BooleanValue(default=False), group='herald')
        self.setattr_argument("enable_force_herald",    BooleanValue(default=False), group='herald')
        self.setattr_argument("force_herald_threshold", NumberValue(default=46, precision=0, step=10, min=0, max=10000), group='herald')

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

        # subsequences
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_rap, ram_addr_start=0, num_samples=500,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # convert deshelving sweep parameters
        freq_deshelve_ftw_list =    np.array([self.repump_qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_deshelve_mhz_list])
        ampl_deshelve_pct_list =    np.array([self.repump_qubit.frequency_to_ftw(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_deshelve_pct_list])
        time_deshelve_mu_list =     np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_deshelve_us_list])

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
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(
            freq_deshelve_ftw_list,
            ampl_deshelve_pct_list,
            time_deshelve_mu_list
        ), -1).reshape(-1, 3)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # configure RAP & record onto DMA
        self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
        with self.core_dma.record('RAP_SUBSEQUENCE'):
            # ensure correct att set (RAP subseq doesn't do this for us)
            self.qubit.set_att_mu(self.att_rap_mu)
            self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        # instantiate relevant variables
        counts_her = -1     # store heralded counts
        # retrieve relevant DMA sequences.handles
        dma_handle_rap = self.core_dma.get_handle('RAP_SUBSEQUENCE')
        self.core.break_realtime()

        # MAIN EXECUTION LOOP
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_deshelve_ftw = config_vals[0]
                ampl_deshelve_asf = config_vals[1]
                time_deshelve_mu =  config_vals[2]
                self.core.break_realtime()

                # set up deshelving/quench/854nm DDS
                self.repump_qubit.set_mu(freq_deshelve_ftw, asf=ampl_deshelve_asf, profile=self.profile_854_target)
                self.core.break_realtime()

                '''PREPARE TARGET SPIN STATE'''
                while True:

                    '''INITIALIZE'''
                    # initialize ion in S-1/2 bright state via Doppler + optical pumping
                    self.initialize_subsequence.run_dma()
                    # run RAP on carrier into D-5/2 dark state
                    self.core_dma.playback_handle(dma_handle_rap)

                    # optional: herald ion via state-dependent fluorescence
                    if self.enable_herald:
                        self.readout_subsequence.run_dma()
                        # make sure to turn beams off after sequence
                        self.pump.off()
                        self.repump_cooling.off()

                        # optional: force heralding (i.e. run until we succeed w/RAP)
                        if self.enable_force_herald:
                            counts_her = self.readout_subsequence.fetch_count()

                            # bright state - try again
                            if counts_her > self.force_herald_threshold:
                                self.core.break_realtime()
                                continue

                            # otherwise, dark state - proceed (and add minor slack)
                            at_mu(self.core.get_rtio_counter_mu() + self.time_force_herald_slack_mu)

                    # force break loop by default
                    break

                '''APPLY DESHELVING'''
                # note: no need to manually pulse IO_UPDATE b/c device class does it for us
                self.repump_qubit.set_profile(self.profile_854_target)
                self.repump_qubit.on()
                delay_mu(time_deshelve_mu)
                self.repump_qubit.off()

                '''READOUT & STORE RESULTS'''
                # read out to verify deshelving
                self.readout_subsequence.run_dma()

                # retrieve heralded measurement ONLY IF THERE EXISTS ONE
                if self.enable_herald and not self.enable_force_herald:
                    counts_her = self.readout_subsequence.fetch_count()

                # retrieve readout results & update dataset
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(freq_deshelve_ftw,
                                    counts_res,
                                    counts_her,
                                    ampl_deshelve_asf,
                                    time_deshelve_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

