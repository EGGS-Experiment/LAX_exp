import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, SidebandCool, RabiFlop, Readout


class SidebandCooling2(LAXExperiment, Experiment):
    """
    Experiment: Sideband Cooling 2

    Measures temperature after a given number of RSB pulses.
    """

    name = 'Sideband Cooling 2'

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                            NumberValue(default=15, ndecimals=0, step=1, min=1, max=10000))

        # sideband cooling readout
        self.setattr_argument("freq_rsb_scan_mhz",                      Scannable(
                                                                            default=CenterScan(103.655, 0.04, 0.001),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group='sideband_readout')
        self.setattr_argument("freq_bsb_scan_mhz",                      Scannable(
                                                                            default=CenterScan(105.271, 0.04, 0.001),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group='sideband_readout')
        self.setattr_argument("time_readout_pipulse_us",                NumberValue(default=250, ndecimals=5, step=1, min=1, max=10000), group='sideband_readout')
        self.setattr_argument("ampl_readout_pipulse_pct",               NumberValue(default=16.1, ndecimals=5, step=1, min=1, max=100), group='sideband_readout')

        # get relevant devices
        self.setattr_device('qubit')

        # get subsequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.sidebandcool_subsequence =                             SidebandCool(self)
        self.rabiflop_subsequence =                                 RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.readout_subsequence =                                  Readout(self)

    def prepare_experiment(self):
        # convert frequencies to machine units
        self.freq_qubit_scan_ftw =                                  np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz])

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(list(self.freq_qubit_scan_mhz)), 2)))
        self.setattr_dataset('results')


    # MAIN SEQUENCE
    @kernel
    def initialize_experiment(self):
        self.core.reset()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.rabiflop_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            self.core.break_realtime()

            # sweep frequency
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=0x1FFF, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # rabi flop
                self.rabiflop_subsequence.run_dma()

                # read out
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

            self.set_dataset('management.completion_pct', (trial_num + 1) / self.repetitions * 100., broadcast=True, persist=True, archive=False)

    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, counts):
        self.results[self._result_iter] = np.array([self.qubit.ftw_to_frequency(freq_ftw) / MHz, counts])
        self._result_iter += 1
