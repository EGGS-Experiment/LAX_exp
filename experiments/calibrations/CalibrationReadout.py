import numpy as np
from artiq.experiment import *
from.artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RescueIon


class CalibrationReadout(LAXExperiment, Experiment):
    """
    Calibration: Readout

    Calibrate/measure readout (i.e. 397nm + 866nm) pulse parameters.
    """
    name = 'Calibration Readout'
    kernel_invariants = {
        # hardware objects
        'initialize_subsequence', 'rescue_subsequence',

        # hardware parameters
        'time_background_pump_mu',

        # experiment values
        'profile_397_readout', 'profile_866_readout', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments/configuration
        self.setattr_argument("repetitions",        NumberValue(default=40, precision=0, step=1, min=1, max=100000))

        # todo: actually implement the doppler thing lol
        self.setattr_argument("enable_doppler",     BooleanValue(default=True), group=None,
                              tooltip="Enable doppler cooling before each data point.")
        self.setattr_argument("time_readout_us_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([10., 20.]),
                                                            RangeScan(10, 100, 91, randomize=True),
                                                        ],
                                                        global_min=1, global_max=2000, global_step=10,
                                                        unit="us", scale=1, precision=3
                                                    ), group=self.name)

        # readout scan parameters - 397nm
        self.setattr_argument("freq_397_readout_mhz_list", Scannable(
                                                        default=[
                                                            CenterScan(110, 3, 0.2, randomize=True),
                                                            ExplicitScan([110.]),
                                                            RangeScan(90., 130., 80, randomize=True),
                                                        ],
                                                        global_min=60, global_max=200, global_step=0.01,
                                                        unit="MHz", scale=1, precision=6
                                                    ), group='397nm')
        self.setattr_argument("ampl_397_readout_pct_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([40., 45.]),
                                                            RangeScan(1, 50, 50, randomize=True),
                                                        ],
                                                        global_min=0.1, global_max=50, global_step=10,
                                                        unit="pct", scale=1, precision=3
                                                    ), group='397nm')

        # readout scan parameters - 866nm
        self.setattr_argument("freq_866_readout_mhz_list", Scannable(
                                                        default=[
                                                            CenterScan(110, 3, 0.2, randomize=True),
                                                            ExplicitScan([110.]),
                                                            RangeScan(90., 130., 80, randomize=True),
                                                        ],
                                                        global_min=60, global_max=200, global_step=0.01,
                                                        unit="MHz", scale=1, precision=6
                                                    ), group='866nm')
        self.setattr_argument("ampl_866_readout_pct_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([21., 25.]),
                                                            RangeScan(1, 50, 50, randomize=True),
                                                        ],
                                                        global_min=0.1, global_max=50, global_step=10,
                                                        unit="pct", scale=1, precision=3
                                                    ), group='866nm')

        '''STATE PREPARATION - CONFIGURATION'''
        # allocate DDS profiles
        self.profile_397_readout =  1
        self.profile_866_readout =  1

        # relevant devices
        self.setattr_device('pmt')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')

        # subsequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # convert deshelving sweep parameters
        time_readout_mu_list =      np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_readout_us_list])
        freq_397_readout_ftw_list = np.array([self.pump.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_397_readout_mhz_list])
        ampl_397_readout_asf_list = np.array([self.pump.amplitude_to_asf(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_397_readout_pct_list])
        freq_866_readout_ftw_list = np.array([self.repump_cooling.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_866_readout_mhz_list])
        ampl_866_readout_asf_list = np.array([self.repump_cooling.amplitude_to_asf(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_866_readout_pct_list])

        # other values
        self.time_background_pump_mu = self.core.seconds_to_mu(50 * us)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(
            time_readout_mu_list,
            freq_397_readout_ftw_list, ampl_397_readout_asf_list,
            freq_866_readout_ftw_list, ampl_866_readout_asf_list
        ), -1).reshape(-1, 5)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                7)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # MAIN EXECUTION LOOP
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                time_readout_mu =       config_vals[0]
                freq_397_readout_ftw =  np.int32(config_vals[1])
                ampl_397_readout_asf =  np.int32(config_vals[2])
                freq_866_readout_ftw =  np.int32(config_vals[3])
                ampl_866_readout_asf =  np.int32(config_vals[4])
                self.core.break_realtime()

                # set up readout (i.e. 397nm and 866nm) DDSs
                self.pump.set_mu(freq_397_readout_ftw, asf=ampl_397_readout_asf, profile=self.profile_397_readout,
                                 phase_mode=PHASE_MODE_CONTINUOUS)
                self.repump_cooling.set_mu(freq_866_readout_ftw, asf=ampl_866_readout_asf, profile=self.profile_866_readout,
                                           phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(10000)


                '''INITIALIZE'''
                # initialize ion in S-1/2 bright state via Doppler + optical pumping
                self.initialize_subsequence.run_dma()

                '''READOUT & STORE RESULTS'''
                # set readout waveform & configure beams
                self.pump.readout()
                self.pump.on()
                self.repump_cooling.on()

                # read out signal
                self.pmt.count(time_readout_mu)

                # ensure ion pumped into D-3/2 state
                self.repump_cooling.off()
                delay_mu(self.time_background_pump_mu)
                # read out background
                self.pmt.count(time_readout_mu)

                # retrieve readout results & update dataset
                counts_sig = self.pmt.fetch_count()
                counts_bgr = self.pmt.fetch_count()
                self.update_results(time_readout_mu,
                                    counts_sig,
                                    counts_bgr,
                                    freq_397_readout_ftw,
                                    ampl_397_readout_asf,
                                    freq_866_readout_ftw,
                                    ampl_866_readout_asf,
                                    )
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

