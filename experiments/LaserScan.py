import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon


class LaserScan(LAXExperiment, Experiment):
    """
    Experiment: Laser Scan

    Gets the number of counts as a function of frequency for a fixed time.
    """
    name = 'Laser Scan'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=50, ndecimals=0, step=1, min=1, max=10000))

        # scan parameters
        self.setattr_argument("freq_qubit_scan_mhz",                Scannable(
                                                                        default=CenterScan(103.345, 0.05, 0.0005, randomize=True),
                                                                        global_min=60, global_max=200, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ), group=self.name)
        self.setattr_argument("time_qubit_us",                      NumberValue(default=5000, ndecimals=5, step=1, min=1, max=10000000), group=self.name)
        self.setattr_argument("att_qubit_db",                       NumberValue(default=28, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)

        # relevant devices
        self.setattr_device('qubit')

        # subsequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.rabiflop_subsequence =                                 RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # convert frequencies to machine units
        self.freq_qubit_scan_ftw =                                  np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz])

        # convert attenuation to machine units
        self.att_qubit_mu =                                         att_to_mu(self.att_qubit_db * dB)


    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_qubit_scan_mhz),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # reduce attenuation/power of qubit beam to resolve lines
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.rabiflop_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            self.core.break_realtime()

            # sweep frequency
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=0x1FFF)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # rabi flop
                self.rabiflop_subsequence.run_dma()

                # read out
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)


    # ANALYSIS
    def analyze(self):
        """
        Fit data and guess potential spectral peaks.
        """
        # todo: separate results readout frequency (bsb vs rsb)
        # todo: calculate count threshold and binarize
        # todo: find peaks
        # todo: save result to dataset
        # todo: print result
        pass
