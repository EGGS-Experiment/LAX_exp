import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout


class QubitAlignment(LAXExperiment, Experiment):
    """
    Utility: Qubit Alignment

    Excite a dark-state qubit transition for a short, fixed time
    to examine the qubit beam alignment.
    """
    name = 'Qubit Alignment'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # general
        self.setattr_argument('time_total_s',           NumberValue(default=500, ndecimals=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('samples_per_point',      NumberValue(default=100, ndecimals=0, step=10, min=1, max=500), group='timing')

        # qubit
        self.setattr_argument('time_qubit_us',          NumberValue(default=3.5, ndecimals=3, step=10, min=0.1, max=100000), group='qubit')
        self.setattr_argument("freq_qubit_mhz",         NumberValue(default=102.9616, ndecimals=5, step=1, min=1, max=10000), group='qubit')
        self.setattr_argument("att_qubit_db",           NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group='qubit')
        self.setattr_argument('counts_threshold',       NumberValue(default=60, ndecimals=0, step=10, min=1, max=1000), group='qubit')


        # instantiate subsequences
        self.initialize_subsequence =                   InitializeQubit(self)
        self.rabiflop_subsequence =                     RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.readout_subsequence =                      Readout(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        # calculate time taken for a single point
        self.time_per_point_us =                ((self.initialize_subsequence.time_repump_qubit_mu
                                                 + self.initialize_subsequence.time_doppler_cooling_mu
                                                 + self.initialize_subsequence.time_spinpol_mu
                                                 + self.readout_subsequence.time_readout_mu) / 1000
                                                 + self.time_qubit_us)

        # get relevant timings and calculate the number of repetitions
        self.repetitions =                      round(self.time_total_s / (self.samples_per_point * self.time_per_point_us * us))

        # declare loop iterators and holder variables ahead of time to reduce overhead
        self._iter_repetitions =                np.arange(self.repetitions)
        self._iter_loop =                       np.arange(self.samples_per_point)
        self._state_array =                     np.zeros(self.samples_per_point, dtype=np.int32)

        # prepare datasets for storing counts
        self.set_dataset('temp.qubit_align._tmp_counts_x',  np.zeros(self.repetitions), broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.qubit_align._tmp_counts_y',  np.zeros(self.repetitions), broadcast=True, persist=False, archive=False)

        # convert qubit parameters
        self.freq_qubit_ftw =                   hz_to_ftw(self.freq_qubit_mhz * MHz)
        self.att_qubit_mu =                     att_to_mu(self.att_qubit_db * dB)

    @property
    def results_shape(self):
        return (self.repetitions, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # prepare qubit beam for readout
        self.qubit.set_profile(0)
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit.set_mu(self.freq_qubit_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)
        self.core.break_realtime()

        # record alignment sequence
        with self.core_dma.record('_QUBIT_ALIGNMENT'):
            # initialize ion in S-1/2 state
            self.initialize_subsequence.run()

            # rabiflop
            self.rabiflop_subsequence.run()

            # read out ion state
            self.readout_subsequence.run()
        
        
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()

        # retrieve DMA handles for qubit alignment
        _handle_alignment = self.core_dma.get_handle('_QUBIT_ALIGNMENT')
        self.core.break_realtime()


        # MAIN LOOP
        for i in self._iter_repetitions:
            self.core.break_realtime()

            # average readings over N samples
            for j in self._iter_loop:

                # run qubit alignment sequence
                self.core_dma.playback_handle(_handle_alignment)

                # get PMT counts
                self._state_array[j] = self.pmt.fetch_count()
                delay_mu(100000)

            # update dataset
            with parallel:
                self.update_results(i, self._state_array)
                self.core.break_realtime()


    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num, state_array):
        # convert raw counts into a probability
        dstate_probability = np.sum((state_array < self.counts_threshold) * 1.) / self.samples_per_point

        # update datasets
        self.mutate_dataset('temp.qubit_align._tmp_counts_x', self._result_iter, iter_num * (self.samples_per_point * self.time_per_point_us * us))
        self.mutate_dataset('temp.qubit_align._tmp_counts_y', self._result_iter, dstate_probability)

        # update completion monitor
        self.set_dataset('management.completion_pct',
                         round(100. * self._result_iter / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += 1

    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
