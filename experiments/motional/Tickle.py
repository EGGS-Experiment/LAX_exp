from artiq.experiment import *
from numpy import array, zeros, int32, int64
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, DopplerRecooling, RescueIon


class Tickle(LAXExperiment, Experiment):
    """
    Experiment: Tickle

    # todo: document
    """
    name = 'Tickle'
    kernel_invariants = {
        # hardware parameters
        "att_tickle_mu", "config_experiment_list", "profile_tickle",

        # subsequences
        'initialize_subsequence', 'readout_counts_subsequence', 'readout_timestamped_subsequence',
        'rescue_subsequence'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # readout configuration
        self.setattr_argument("readout_type",   EnumerationValue(['Counts', 'Timestamped'], default='Timestamped'), group=self.name)

        # tickle configuration
        self.setattr_argument("tickle_source",  EnumerationValue(['Parametric', 'Dipole'], default='Dipole'), group=self.name)
        self.setattr_argument("freq_tickle_khz_list",   Scannable(
                                                            default=[
                                                                CenterScan(1252.6, 20, 0.1, randomize=True),
                                                                ExplicitScan([1500]),
                                                            ],
                                                            global_min=10, global_max=400000, global_step=100,
                                                            unit="kHz", scale=1, precision=3
                                                        ), group=self.name)
        self.setattr_argument("ampl_tickle_pct_list",   Scannable(
                                                            default=[
                                                                ExplicitScan([35.]),
                                                                RangeScan(0, 100, 10, randomize=True),
                                                            ],
                                                            global_min=0.1, global_max=100.0, global_step=10,
                                                            unit="pct", scale=1, precision=2
                                                        ), group=self.name)
        self.setattr_argument("time_tickle_us_list",    Scannable(
                                                            default=[
                                                                ExplicitScan([1000.]),
                                                                RangeScan(1000, 2000, 11, randomize=True),
                                                            ],
                                                            global_min=100, global_max=1000000, global_step=100,
                                                            unit="us", scale=1, precision=0
                                                        ), group=self.name)
        self.setattr_argument("att_tickle_db",  NumberValue(default=10, precision=1, step=0.5, min=0., max=31.5), group=self.name)

        # get necessary devices
        self.setattr_device('dds_parametric')
        self.setattr_device('dds_dipole')

        # assign DDS profiles
        self.profile_tickle = 0

        # subsequences
        self.initialize_subsequence =           InitializeQubit(self)
        self.readout_counts_subsequence =       Readout(self)
        self.readout_timestamped_subsequence =  DopplerRecooling(self)
        self.rescue_subsequence =               RescueIon(self)

    def prepare_experiment(self):
        # select desired readout method
        if self.readout_type == 'Counts':           self.readout_subsequence = self.readout_counts_subsequence
        elif self.readout_type == 'Timestamped':    self.readout_subsequence = self.readout_timestamped_subsequence

        # select desired tickle subsequence based on input arguments
        if self.tickle_source == 'Parametric':  self.dds_tickle = self.dds_parametric
        elif self.tickle_source == 'Dipole':    self.dds_tickle = self.dds_dipole

        # convert tickle parameters to machine units
        self.att_tickle_mu =        att_to_mu(self.att_tickle_db * dB)
        freq_tickle_ftw_list =  [self.dds_tickle.frequency_to_ftw(freq_khz * kHz)
                                 for freq_khz in self.freq_tickle_khz_list]
        ampl_tickle_asf_list =  [self.dds_tickle.amplitude_to_asf(ampl_pct / 100.)
                                 for ampl_pct in self.ampl_tickle_pct_list]
        time_tickle_mu_list =   [self.core.seconds_to_mu(time_us * us)
                                 for time_us in self.time_tickle_us_list]

        # create an array of values for the experiment to sweep
        # (i.e. tickle frequency, tickle time, readout FTW)
        self.config_experiment_list = create_experiment_config(
            freq_tickle_ftw_list, ampl_tickle_asf_list, time_tickle_mu_list,
            shuffle_config=True, config_type=int64
        )

        # configure timestamped counts for doppler recooling
        if self.readout_type == 'Timestamped':
            self.readout_subsequence = self.readout_timestamped_subsequence
            self.set_dataset('timestamped_counts', zeros((self.repetitions * len(self.config_experiment_list),
                                                          self.readout_subsequence.num_counts), dtype=int64))
            self.setattr_dataset('timestamped_counts')

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)


    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # ensure DMA sequences use profile 0
        self.dds_tickle.set_profile(self.profile_tickle)
        self.dds_tickle.set_att_mu(self.att_tickle_mu)
        delay_mu(10000)

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_tickle_ftw =   int32(config_vals[0])
                ampl_tickle_asf =   int32(config_vals[1])
                time_tickle_mu =    config_vals[2]

                # configure tickle and qubit readout
                self.core.break_realtime()
                self.dds_tickle.set_mu(freq_tickle_ftw, asf=ampl_tickle_asf,
                                       profile=self.profile_tickle,
                                       phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(8000)

                '''INITIALIZE ION & EXCITE'''
                # initialize ion
                self.initialize_subsequence.run_dma()

                # tickle
                self.dds_tickle.on()
                delay_mu(time_tickle_mu)
                self.dds_tickle.off()

                '''READ OUT AND RECORD RESULTS'''
                # read out results and clean up loop
                self.readout_subsequence.run()
                self.rescue_subsequence.resuscitate()
                self.update_results(
                    freq_tickle_ftw,
                    self.readout_subsequence.fetch_count(),
                    ampl_tickle_asf,
                    time_tickle_mu
                )

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()

    @rpc(flags={"async"})
    def update_results(self, *args) -> TNone:
        # tmp remove
        counts_tmp =    0
        res_tmp =       array([])
        if self.readout_type == "Timestamped":
            counts_tmp = args[1][-1]
            self.mutate_dataset('timestamped_counts', self._result_iter, args[1])
            res_tmp = array([args[0], counts_tmp, args[2], args[3]])
        else:
            res_tmp = array([args])
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

