import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout

# tmp remove
from collections import deque
# tmp remove


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
        self.setattr_argument('time_total_s',           NumberValue(default=500, precision=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('samples_per_point',      NumberValue(default=50, precision=0, step=10, min=15, max=500), group='timing')

        # qubit
        self.setattr_argument('time_qubit_us',          NumberValue(default=3.5, precision=3, step=10, min=0.1, max=100000), group='qubit')
        self.setattr_argument("freq_qubit_mhz",         NumberValue(default=101.4518, precision=5, step=1, min=1, max=10000), group='qubit')
        self.setattr_argument("att_qubit_db",           NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5), group='qubit')
        self.setattr_argument('counts_threshold',       NumberValue(default=46, precision=0, step=10, min=1, max=1000), group='qubit')


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

        # convert qubit parameters
        self.freq_qubit_ftw =                   hz_to_ftw(self.freq_qubit_mhz * MHz)
        self.att_qubit_mu =                     att_to_mu(self.att_qubit_db * dB)

        # tmp remove
        # todo: set the filter up properly
        self._th_wind = 100
        self._th0 = deque(maxlen=self._th_wind)
        # tmp remove

    @rpc
    def initialize_plotting(self) -> TNone:
        """
        Configure datasets and applets for plotting.
        """
        # prepare datasets for storing counts
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        counts_x_arr = np.zeros(self.repetitions) * np.nan
        counts_x_arr[0] = 0
        counts_y_arr = np.zeros(self.repetitions) * np.nan
        counts_y_arr[0] = 0

        # prepare datasets for storing counts
        self.set_dataset('temp.qubit_align.counts_x', counts_x_arr, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.qubit_align.counts_y', counts_y_arr, broadcast=True, persist=False, archive=False)

        # initialize plotting applet
        self.ccb.issue("create_applet", "qubit_alignment",
                       '${artiq_applet}plot_xy temp.qubit_align.counts_y'
                       ' --x temp.qubit_align.counts_x --title "Qubit Alignment"')

    @property
    def results_shape(self):
        return (self.repetitions, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # prepare plotting
        # note: do it here instead of prepare to prevent overriding other experiments
        self.initialize_plotting()
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
        for num_rep in self._iter_repetitions:
            self.core.break_realtime()

            # average readings over N samples
            for num_count in self._iter_loop:

                # run qubit alignment sequence
                self.core_dma.playback_handle(_handle_alignment)

                # get PMT counts
                self._state_array[num_count] = self.pmt.fetch_count()
                delay_mu(10000)

            # update dataset
            self.update_results(num_rep, self._state_array)
            self.core.break_realtime()

            # periodically check termination
            if (num_rep % 10) == 0:
                self.check_termination()
                self.core.break_realtime()


    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num, state_array):
        # convert raw counts into a probability
        # dstate_probability = np.sum((state_array < self.counts_threshold) * 1.) / self.samples_per_point
        # tmp remove
        for val in state_array:
            self._th0.append((val < self.counts_threshold) * 1.)
        # tmp remove


        # update datasets
        self.mutate_dataset('temp.qubit_align.counts_x', self._result_iter, iter_num * (self.samples_per_point * self.time_per_point_us * us))
        # self.mutate_dataset('temp.qubit_align.counts_y', self._result_iter, dstate_probability)
        self.mutate_dataset('temp.qubit_align.counts_y', self._result_iter, np.mean(self._th0))

        # update completion monitor
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._result_iter / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += 1

