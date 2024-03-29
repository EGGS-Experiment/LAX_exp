import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.SidebandCooling as SidebandCooling
from LAX_exp.system.subsequences import TickleDDS, TickleFastDDS, TickleFastPhaser


class QLMSRabiRDX(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Rabi - RDX

    Quantum Logic Mass Spectroscopy - Rabi Flopping (RDX)
    Cool the ions to the ground state of motion via sideband cooling,
    then apply a tickle to create a coherent state to be read out via RSB/BSB comparison.
    *** todo: add documentation about being fast and time sweeping
    """
    name = 'QLMSRabiRDX'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # QLMS configuration
        self.setattr_argument("freq_qlms_rabi_mhz_list",                            Scannable(
                                                                                        default=[
                                                                                            ExplicitScan([82]),
                                                                                            CenterScan(82, 2, 0.5, randomize=True)
                                                                                        ],
                                                                                        global_min=0, global_max=400, global_step=1,
                                                                                        unit="MHz", scale=1, ndecimals=5
                                                                                    ), group=self.name)
        self.setattr_argument("phase_qlms_rabi_turns_list",                         Scannable(
                                                                                        default=[
                                                                                            ExplicitScan([0.]),
                                                                                            RangeScan(0, 1.0, 6, randomize=True)
                                                                                        ],
                                                                                        global_min=0.0, global_max=1.0, global_step=1,
                                                                                        unit="turns", scale=1, ndecimals=3
                                                                                    ), group=self.name)
        self.setattr_argument("time_qlms_rabi_ns_list",                             Scannable(
                                                                                        default=[
                                                                                            ExplicitScan([8]),
                                                                                            RangeScan(8, 250, 38, randomize=True)
                                                                                        ],
                                                                                        global_min=40, global_max=10000000, global_step=40,
                                                                                        unit="ns", scale=1, ndecimals=0
                                                                                    ), group=self.name)

        # set up tickle source selection
        self.setattr_argument("tickle_source",                                      EnumerationValue(['DDS', 'Phaser'], default='DDS'), group='ticklefast')
        self.tickle_subsequence_dds =                                               TickleFastDDS(self)
        self.tickle_subsequence_phaser =                                            TickleFastPhaser(self)

        # get necessary devices
        self.setattr_device('dds_parametric')

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list =                                   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list

        # select desired tickle subsequence based on input arguments
        if self.tickle_source == 'DDS':                                         self.tickle_subsequence = self.tickle_subsequence_dds
        elif self.tickle_source == 'Phaser':                                    self.tickle_subsequence = self.tickle_subsequence_phaser

        # convert QLMS modulation to machine units
        self.freq_qlms_rabi_ftw_list =                                          np.array([hz_to_ftw(freq_mhz * MHz)
                                                                                          for freq_mhz in self.freq_qlms_rabi_mhz_list])
        self.phase_qlms_rabi_pow_list =                                         np.array([self.dds_parametric.turns_to_pow(phase_turns)
                                                                                          for phase_turns in self.phase_qlms_rabi_turns_list])
        self.time_qlms_rabi_mu_list =                                           np.array([self.core.seconds_to_mu(time_ns * ns)
                                                                                          for time_ns in self.time_qlms_rabi_ns_list])


        # ROUND (not truncate) values to nearest switching time multiple
        self.time_qlms_rabi_mu_list = np.array(list(self.time_qlms_rabi_ns_list))
        if self.tickle_source == 'DDS':
            self.time_qlms_rabi_mu_list = np.int64(np.round(self.time_qlms_rabi_mu_list + 7)) & ~0x7
        elif self.tickle_source == 'Phaser':
            self.time_qlms_rabi_mu_list = np.array([40 * max(1.0, time_mu)
                                                    for time_mu in np.round(self.time_qlms_rabi_mu_list / 40)],
                                                   dtype = np.int64)
            # self.time_qlms_rabi_mu_list = np.int64(np.round(np.floor(self.time_qlms_rabi_mu_list / 40 + 0.501)) * 40)

        # remove duplicate elements to stop us from wasting time
        self.time_qlms_rabi_mu_list = np.unique(self.time_qlms_rabi_mu_list)
        print('\ttimes: {}'.format(self.time_qlms_rabi_mu_list))
        print('\t\tlen: {}'.format(len(self.time_qlms_rabi_mu_list)))


        # create an array of values for the experiment to sweep
        # (i.e. tickle frequency, tickle time, readout FTW)
        self.config_qlms_rabi_list =                                            np.stack(np.meshgrid(self.freq_qlms_rabi_ftw_list,
                                                                                                     self.phase_qlms_rabi_pow_list,
                                                                                                     self.time_qlms_rabi_mu_list,
                                                                                                     self.freq_sideband_readout_ftw_list),
                                                                                         -1).reshape(-1, 4)
        np.random.shuffle(self.config_qlms_rabi_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_qlms_rabi_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):
            for config_vals in self.config_qlms_rabi_list:

                # extract values from config list
                freq_qlms_ftw =     np.int32(config_vals[0])
                phase_qlms_pow =    np.int32(config_vals[1])
                time_qlms_mu =      config_vals[2]
                freq_readout_ftw =  np.int32(config_vals[3])
                self.core.break_realtime()

                # configure tickle and qubit readout
                self.tickle_subsequence.configure(freq_qlms_ftw, phase_qlms_pow, time_qlms_mu)
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # QLMS fast tickle
                self.tickle_subsequence.run()

                # sideband readout
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(),
                                        freq_qlms_ftw, phase_qlms_pow, time_qlms_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

    def analyze(self):
        pass
