import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon

_THRESHOLD1 = 56
_THRESHOLD2 = 282


class RabiFloppingDynamic(LAXExperiment, Experiment):
    """
    Experiment: Rabi Flopping Dynamic

    Measures ion fluorescence vs 729nm pulse time and frequency. dynamic
    """
    name = 'Rabi Flopping Dynamic'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=30, ndecimals=0, step=1, min=1, max=10000))

        # rabi flopping arguments
        self.setattr_argument("time_rabi_us_list",                  Scannable(
                                                                        default=RangeScan(0, 50, 101, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ), group=self.name)
        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=103.855, ndecimals=5, step=1, min=1, max=10000), group=self.name)
        self.setattr_argument("att_qubit_db",                       NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)


        # get devices
        self.setattr_device('qubit')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # convert rabi flopping arguments
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [seconds_to_mu((max_time_us - time_us) * us), seconds_to_mu(time_us * us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])
        self.freq_rabiflop_ftw =                                    hz_to_ftw(self.freq_rabiflop_mhz * MHz)

        # convert attenuation to machine units
        self.att_qubit_mu =                                         att_to_mu(self.att_qubit_db * dB)

        # tmp remove
        self.datatmp = np.zeros(len(self.time_rabiflop_mu_list))
        self.set_dataset('datatmp_holder', np.zeros((len(self.time_rabiflop_mu_list), 2)))
        self.setattr_dataset('datatmp_holder')
        self._iter_tmp = 0
        # tmp remove


    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # set qubit beam parameters
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # sweep time
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # wait given time
                delay_mu(time_rabi_pair_mu[0])

                # rabi flopping w/qubit laser
                self.qubit.on()
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.off()

                # do readout
                self.readout_subsequence.run_dma()

                # tmp remove
                counts = self.readout_subsequence.fetch_count()
                # tmp remove

                # update dataset
                with parallel:
                    self.update_results(time_rabi_pair_mu[1], counts)
                    self.update_results2(time_rabi_pair_mu[1], counts)
                    self.core.break_realtime()


            # tmp remove
            with parallel:
                self.update_results_dynamic(trial_num)
                self.core.break_realtime()
            # tmp remove

    # tmp remove
    @rpc(flags={"async"})
    def update_results2(self, time, counts):
        self.mutate_dataset('datatmp_holder', self._iter_tmp, np.array([time, counts]))
        self._iter_tmp = (self._iter_tmp + 1) % len(self.time_rabiflop_mu_list)

    @rpc(flags={"async"})
    def update_results_dynamic(self, n):
        """
        todo: document
        """
        tmp = np.array(self.datatmp_holder)

        # get counts
        counts = tmp[np.argsort(tmp, axis=0)[:, 0]][:, 1]
        probability = (counts >= _THRESHOLD1) * 0.5 + (counts >= _THRESHOLD2) * 0.5

        # calculate moving average and update
        # self.datatmp = np.array((self.datatmp * (n) + probability) / (n + 1))
        alpha = 0.25
        self.datatmp = np.array(alpha * probability + (1 - alpha) * self.datatmp)
        self.set_dataset('datatmp', self.datatmp, broadcast=True,persist=True,archive=True)
    # tmp remove
