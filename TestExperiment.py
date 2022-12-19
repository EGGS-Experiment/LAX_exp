import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment
from LAX_exp.utilities.conversions import *
from LAX_exp.system.subsequences import RabiFlop, DopplerCool, Readout


class TestExperiment(LAXExperiment, Experiment):
    """
    Test Experiment
    Used for testing. Runs a Laser Scan.
    """

    name = 'Laser Scan 2'

    def build_experiment(self):
        # timing
        self.setattr_argument("time_729_us",                    NumberValue(default=1, ndecimals=5, step=1, min=1, max=10000000))

        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(
                                                                    default=RangeScan(104.24, 104.96, 2, randomize=True),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ))

        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(
                                                                    default=LinearScan(104.24, 104.96, 2, randomize=True),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ))

    def prepare_experiment(self):
        # get 729 beam
        self.setattr_device('qubit')
        self.setattr_device('pmt')

        # prepare sequences
        self.cooling_subsequence =                              DopplerCool(self)
        self.rabiflop_subsequence =                             RabiFlop(self, time_rabiflop_us=self.time_729_us)
        self.readout_subsequence =                              Readout(self)

        # dataset
        self.set_dataset('results',                             np.zeros((self.repetitions * len(list(self.freq_qubit_scan_mhz)), 2)))
        self.setattr_dataset('results')

        # prepare frequencies for sweep
        self.freq_qubit_scan_ftw =                              [mhz_to_ftw(freq_mhz) for freq_mhz in self.freq_qubit_scan_mhz]
        self.pmt = self.get_device("pmt")


    # PREPARE MAIN SEQUENCE
    def run_initialize(self):
        # record pulse sequence onto DMA for speed & accuracy
        dma_handle = self._record_dma()
        setattr(self, 'dma_handle', dma_handle)

    @kernel(flags='fast-math')
    def _record_dma(self):
        # record DMA sequence
        with self.core_dma.record('main_sequence'):
            self.cooling_subsequence.run()
            self.rabiflop_subsequence.run()
            self.readout_subsequence.run()
        self.core.break_realtime()

        # get and return DMA handle
        handle = self.core_dma.get_handle('main_sequence')
        self.core.break_realtime()
        return handle


    # MAIN LOOP
    @kernel
    def run_main(self):
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep frequency
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=0x1FFF)
                self.core.break_realtime()

                # run main sequence
                self.core_dma.playback_handle(self.dma_handle)

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, self.pmt.fetch_count())
                    self.core.break_realtime()
