import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, DopplerRecooling, RescueIon


class Tickle(LAXExperiment, Experiment):
    """
    Experiment: Tickle

    # todo: document
    """
    name = 'Tickle'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                    NumberValue(default=50, ndecimals=0, step=1, min=1, max=10000))

        # readout configuration
        self.setattr_argument("readout_type",                   EnumerationValue(['Counts', 'Timestamped'], default='Timestamped'), group=self.name)

        # tickle configuration
        self.setattr_argument("tickle_source",                  EnumerationValue(['Parametric', 'Dipole'], default='Dipole'), group=self.name)
        self.setattr_argument("freq_tickle_khz_list",           Scannable(
                                                                        default=[
                                                                            CenterScan(1252.6, 20, 0.1, randomize=True),
                                                                            ExplicitScan([1500]),
                                                                        ],
                                                                        global_min=10, global_max=400000, global_step=100,
                                                                        unit="kHz", scale=1, ndecimals=3
                                                                    ), group=self.name)
        self.setattr_argument("ampl_tickle_pct_list",           Scannable(
                                                                        default=[
                                                                            ExplicitScan([35.]),
                                                                            RangeScan(0, 100, 10, randomize=True),
                                                                        ],
                                                                        global_min=0.1, global_max=100.0, global_step=10,
                                                                        unit="pct", scale=1, ndecimals=2
                                                                    ), group=self.name)
        self.setattr_argument("time_tickle_us_list",            Scannable(
                                                                        default=[
                                                                            ExplicitScan([1000.]),
                                                                            RangeScan(1000, 2000, 11, randomize=True),
                                                                        ],
                                                                        global_min=100, global_max=1000000, global_step=100,
                                                                        unit="us", scale=1, ndecimals=0
                                                                    ), group=self.name)
        self.setattr_argument("att_tickle_db",                  NumberValue(default=10, ndecimals=1, step=0.5, min=2, max=31.5), group=self.name)

        # get necessary devices
        self.setattr_device('dds_parametric')
        self.setattr_device('dds_dipole')

        # subsequences
        self.initialize_subsequence =               InitializeQubit(self)
        self.readout_counts_subsequence =           Readout(self)
        self.readout_timestamped_subsequence =      DopplerRecooling(self)
        self.rescue_subsequence =                   RescueIon(self)

    def prepare_experiment(self):
        # select desired readout method
        if self.readout_type == 'Counts':           self.readout_subsequence = self.readout_counts_subsequence
        elif self.readout_type == 'Timestamped':    self.readout_subsequence = self.readout_timestamped_subsequence

        # select desired tickle subsequence based on input arguments
        if self.tickle_source == 'Parametric':      self.dds_tickle = self.dds_parametric
        elif self.tickle_source == 'Dipole':        self.dds_tickle = self.dds_dipole

        # convert tickle parameters to machine units
        self.att_tickle_mu =                        att_to_mu(self.att_tickle_db * dB)
        self.freq_tickle_ftw_list =                 np.array([hz_to_ftw(freq_khz * kHz) for freq_khz in self.freq_tickle_khz_list])
        self.ampl_tickle_asf_list =                 np.array([pct_to_asf(ampl_pct) for ampl_pct in self.ampl_tickle_pct_list])
        self.time_tickle_mu_list =                  np.array([self.core.seconds_to_mu(time_us * us) for time_us in self.time_tickle_us_list])

        # create an array of values for the experiment to sweep
        # (i.e. tickle frequency, tickle time, readout FTW)
        self.config_tickle_list =                   np.stack(np.meshgrid(self.freq_tickle_ftw_list,
                                                                         self.ampl_tickle_asf_list,
                                                                         self.time_tickle_mu_list),
                                                             -1).reshape(-1, 3)
        np.random.shuffle(self.config_tickle_list)

        # tmp remove
        if self.readout_type == 'Timestamped':
            self.readout_subsequence = self.readout_timestamped_subsequence
            self.set_dataset('timestamped_counts', np.zeros((self.repetitions * len(self.config_tickle_list), self.readout_subsequence.num_counts), dtype=np.int64))
            self.setattr_dataset('timestamped_counts')
        # tmp remove

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_tickle_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # ensure DMA sequences use profile 0
        self.dds_tickle.set_profile(0)
        self.dds_tickle.set_att_mu(self.att_tickle_mu)

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):
            for config_vals in self.config_tickle_list:

                # extract values from config list
                freq_tickle_ftw =   np.int32(config_vals[0])
                ampl_tickle_asf =   np.int32(config_vals[1])
                time_tickle_mu =    config_vals[2]
                self.core.break_realtime()

                # configure tickle and qubit readout
                self.dds_tickle.set_mu(freq_tickle_ftw, asf=ampl_tickle_asf, profile=0)
                self.core.break_realtime()

                # initialize ion
                self.initialize_subsequence.run_dma()

                # tickle
                self.dds_tickle.on()
                delay_mu(time_tickle_mu)
                self.dds_tickle.off()

                # fluorescence readout
                self.readout_subsequence.run()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_tickle_ftw,
                        self.readout_subsequence.fetch_count(),
                        ampl_tickle_asf,
                        time_tickle_mu
                    )
                    self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        # todo: cleanup

    @rpc(flags={"async"})
    def update_results(self, *args):
        # tmp remove
        res_tmp = np.array([])
        counts_tmp = 0

        if self.readout_type == "Timestamped":
            counts_tmp = args[1][-1]
            self.mutate_dataset('timestamped_counts', self._result_iter, args[1])
            res_tmp = np.array([args[0], counts_tmp, args[2], args[3]])
        else:
            res_tmp = np.array([args])
            counts_tmp = args[1]
        # tmp remove


        # store results in main dataset
        self.mutate_dataset('results', self._result_iter, res_tmp)

        # do intermediate processing
        if (self._result_iter % self._dynamic_reduction_factor) == 0:
            # plot counts in real-time to monitor ion death
            self.mutate_dataset('temp.counts.trace', self._counts_iter, counts_tmp)
            self._counts_iter += 1

            # monitor completion status
            self.set_dataset('management.dynamic.completion_pct', round(self._result_iter * self._completion_iter_to_pct, 3),
                             broadcast=True, persist=True, archive=False)

        # increment result iterator
        self._result_iter += 1

    def analyze(self):
        pass
