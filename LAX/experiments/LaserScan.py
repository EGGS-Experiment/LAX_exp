import numpy as np
from artiq.experiment import *

from LAX_exp.utilities.conversions import *
from LAX_exp.LAX.base_classes import LAXExperiment
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
                                                                    default=RangeScan(104.24, 104.96, 100, randomize=True),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ))

    def prepare_experiment(self):
        # get 729 beam
        self.setattr_device('qubit')
        self.setattr_device('pmt')

        # sequences
        self.cooling_subsequence =                              DopplerCool(self)
        self.rabiflop_subsequence =                             RabiFlop(self, time_rabiflop_us=self.time_729_us)
        self.readout_subsequence =                              Readout(self)

        # dataset
        self.set_dataset('results',                             np.zeros((self.repetitions * len(list(self.freq_qubit_scan_mhz)), 2)))
        self.setattr_dataset('results')

        # tmp remove
        self.freq_qubit_scan_ftw =                              [mhz_to_ftw(freq_mhz) for freq_mhz in self.freq_qubit_scan_mhz]
        self.pmt = self.get_device("pmt")
        self.pump = self.get_device('pump')
        self.qubit = self.get_device('qubit')

    #@kernel
    def run(self):
        self.core.reset()

        handle = self.thkimde()
        setattr(self, 'handle', handle)

        for trial_num in range(self.repetitions):
            self._run_main()

    @kernel
    def _run_main(self):
        self.core.reset()
        for freq_ftw in self.freq_qubit_scan_ftw:
            # set qubit frequency
            self.qubit.set_mu(freq_ftw, asf=0x1FFF)
            self.core.break_realtime()

            self.core_dma.playback_handle(self.handle)

            # update dataset
            with parallel:
                self.update_dataset(freq_ftw, self.pmt.fetch_count())
                self.core.break_realtime()

    # def run_loop(self):
    #     """
    #     Run the experimental sequence.
    #     """
    #     for freq_mhz in self.freq_qubit_scan_mhz:
    #
    #         # convert value to ftw
    #         freq_ftw = mhz_to_ftw(freq_mhz)
    #
    #         # get and store results
    #         try:
    #             counts = self._run_loop_kernel(freq_ftw)
    #             self.update_dataset(freq_mhz, counts)
    #         except RTIOUnderflow:
    #             self.core.break_realtime()
    #
    # @kernel(flags='fast-math')
    # def _run_loop_kernel(self, freq_ftw):
    #     self.core.break_realtime()
    #
    #     # set qubit frequency
    #     self.qubit.set_mu(freq_ftw, asf=self.qubit.ampl_qubit_asf)
    #     self.core.break_realtime()
    #
    #     # doppler cool
    #     self.cooling_subsequence.run_dma()
    #
    #     # rabi flop
    #     self.rabiflop_subsequence.run_dma()
    #
    #     # readout
    #     self.readout_subsequence.run_dma()
    #     self.core.break_realtime()
    #
    #     # return data
    #     return self.pmt.fetch_count()

    @kernel(flags='fast-math')
    def thkimde(self):
        with self.core_dma.record('thkim'):
            self.cooling_subsequence.run()

            self.rabiflop_subsequence.run()

            self.readout_subsequence.run()

        self.core.break_realtime()
        handle = self.core_dma.get_handle('thkim')
        self.core.break_realtime()
        return handle
