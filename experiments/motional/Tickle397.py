from artiq.experiment import *
from numpy import array, zeros, int32, int64

from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, Readout, DopplerRecooling, RescueIon

import os
import labrad
from socket import gethostname


class Tickle397(LAXExperiment, Experiment):
    """
    Experiment: Tickle 397

    # todo: document
    """
    name = 'Tickle 397'
    kernel_invariants = {
        # hardware parameters
        "att_tickle_mu", "config_experiment_list",

        # labrad
        "cxn_wm", "wm",

        # subsequences
        'initialize_subsequence', 'readout_counts_subsequence', 'readout_timestamped_subsequence',
        'rescue_subsequence',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=50, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("wavemeter_chan", NumberValue(default=4, precision=0, step=1, min=1, max=16))
        self.setattr_argument("pzt_397_voltage_v",  Scannable(
                                                            default=[
                                                                ExplicitScan([71.01]),
                                                                RangeScan(68, 71, 100, randomize=True),
                                                            ],
                                                            global_min=10, global_max=100, global_step=1,
                                                            unit="V", scale=1, precision=3
                                                        ))


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
        self.att_tickle_mu =    att_to_mu(self.att_tickle_db * dB)
        freq_tickle_ftw_list =  [self.dds_tickle.frequency_to_ftw(freq_khz * kHz)
                                 for freq_khz in self.freq_tickle_khz_list]
        ampl_tickle_asf_list =  [self.dds_tickle.amplitude_to_asf(ampl_pct / 100.)
                                 for ampl_pct in self.ampl_tickle_pct_list]
        time_tickle_mu_list =   [self.core.seconds_to_mu(time_us * us)
                                 for time_us in self.time_tickle_us_list]

        # check input for serious booboos
        if any(not (20 <= voltage_v <= 90) for voltage_v in self.pzt_397_voltage_v):
            raise ValueError("Invalid fucking 397 pzt voltage you idiot - don't fucking damage it")

        # create an array of values for the experiment to sweep
        self.pzt_397_voltage_v = array(list(self.pzt_397_voltage_v))
        self.config_experiment_list = create_experiment_config(
            freq_tickle_ftw_list, ampl_tickle_asf_list, time_tickle_mu_list,
            shuffle_config=True, config_type=int64
        )

        # tmp remove
        if self.readout_type == 'Timestamped':
            self.readout_subsequence = self.readout_timestamped_subsequence
            self.set_dataset('timestamped_counts', zeros((self.repetitions * len(self.config_experiment_list),
                                                             self.readout_subsequence.num_counts), dtype=int64))
            self.setattr_dataset('timestamped_counts')
        # tmp remove

        # connect to labrad - EGGS
        self.cxn = labrad.connect(os.environ['LABRADHOST'], port=7682, tls_mode='off',
                                  username='', password='lab')
        self.toptica = self.cxn.toptica_server

        # create synchronous connection to wavemeter labrad
        from EGGS_labrad.config.multiplexerclient_config import multiplexer_config
        self.cxn_wm = labrad.connect(multiplexer_config.ip,
                                     name="{:s}_({:s})".format("ARTIQ_EXP", gethostname()),
                                     username="", password=os.environ['LABRADPASSWORD'])
        self.wm = self.cxn_wm.multiplexerserver

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list) * len(self.pzt_397_voltage_v),
                6)


    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # ensure DMA sequences use profile 0
        self.dds_tickle.set_profile(0)
        self.dds_tickle.set_att_mu(self.att_tickle_mu)
        delay_mu(10000)

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for voltage_v in self.pzt_397_voltage_v:
            # update toptica voltage
            voltage_v2 = self._update_toptica(voltage_v)
            self.core.break_realtime()
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()
            delay_mu(10000000)

            for trial_num in range(self.repetitions):
                # get wavemeter frequency - only once per repetition
                wm_freq = self._record_wavemeter()

                for config_vals in self.config_experiment_list:

                    '''CONFIGURE'''
                    # extract values from config list
                    freq_tickle_ftw =   int32(config_vals[0])
                    ampl_tickle_asf =   int32(config_vals[1])
                    time_tickle_mu =    config_vals[2]

                    # configure tickle and qubit readout
                    self.core.break_realtime()
                    self.dds_tickle.set_mu(freq_tickle_ftw, asf=ampl_tickle_asf, profile=0)
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

                    # update dataset
                    self.update_results(
                        freq_tickle_ftw,
                        self.readout_subsequence.fetch_count(),
                        ampl_tickle_asf,
                        time_tickle_mu,
                        voltage_v2,
                        wm_freq
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
            self.set_dataset('management.dynamic.completion_pct',
                             round(self._result_iter * self._completion_iter_to_pct, 3),
                             broadcast=True, persist=True, archive=False)

        # increment result iterator
        self._result_iter += 1


    """
    HELPERS2
    """
    @rpc
    def _update_toptica(self, voltage_v: TFloat) -> TFloat:
        try:
            voltage_v = self.toptica.piezo_set(7, voltage_v)
            # print(voltage_v)
        except Exception as e:
            voltage_v = -1

        return voltage_v

    @rpc
    def _record_wavemeter(self) -> TFloat:
        """
        Attempt to retrieve wavemeter freq for target channel.
        """
        try:
            freq_thz = self.wm.get_frequency(self.wavemeter_chan)
        except Exception as e:
            freq_thz = -1

        return freq_thz

    def analyze(self):
        pass

