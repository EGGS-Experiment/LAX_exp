import numpy as np
from artiq.experiment import *

from LAX_exp.LAX.base_classes import LAXExperiment
from LAX_exp.LAX.devices import
from LAX_exp.LAX.subsequences import RabiFlop, DopplerCool, Readout


class LaserScan2(LAXExperiment):
    """
    729nm Laser Scan2
    Gets the number of counts as a function of frequency for a fixed time.
    """

    name = 'Laser Scan 2'

    def build_experiment(self):
        # timing
        self.setattr_argument("time_729_us",                    NumberValue(default=400, ndecimals=5, step=1, min=1, max=10000000))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(
                                                                    default=RangeScan(104.24, 104.96, 801, randomize=True),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ))

        # get 729 beam
        self.setattr_device('qubit')
        self.setattr_device('pmt')


    def prepare_experiment(self):
        # sequences
        self.cooling_subsequence =                              DopplerCool(self)
        self.rabiflop_subsequence =                             RabiFlop(self, time_rabiflop_us=self.time_729_us)
        self.readout_subsequence =                              Readout(self)

        # dataset
        self.results =                                          np.zeros((self.repetitions * len(list(self.freq_qubit_scan_mhz)), 2))
        self._result_iter =                                     0


    def run_loop(self):
        """
        Run the experimental sequence.
        """
        for freq_mhz in self.freq_qubit_scan_mhz:
            self._run_loop_kernel(freq_mhz)

    @kernel(flags='fast-math')
    def _run_loop_kernel(self, freq_mhz):
        # set qubit frequency
        self.qubit.set(freq_mhz * 1e6, asf=0.5)

        # doppler cool
        self.cooling_subsequence.run_dma()

        # rabi flop
        self.rabiflop_subsequence.run_dma()

        # readout
        self.readout_subsequence.run_dma()

        self.core.break_realtime()
        self.update_dataset()

    @rpc(flags='async')
    def update_dataset(self, freq_mhz, pmt_counts):
        self.results[self.iter] = np.array([freq_mhz, pmt_counts])
        self.iter += 1
