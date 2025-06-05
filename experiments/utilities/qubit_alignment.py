import numpy as np
from artiq.experiment import *

from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, ReadoutAdaptive
# todo: dynamically filter results


class QubitAlignment(LAXExperiment, Experiment):
    """
    Utility: Qubit Alignment

    Excite a dark-state qubit transition for a short, fixed time
    to examine the qubit beam alignment.
    """
    name = 'Qubit Alignment'
    kernel_invariants = {
        # conversions
        "time_per_point_us", "repetitions",
        # hardware parameters etc.
        "freq_qubit_ftw", "att_qubit_mu", "_iter_repetitions", "_iter_loop",
        # subsequences
        "initialize_subsequence", "rabiflop_subsequence", "readout_subsequence"
    }


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # general
        self.setattr_argument('time_total_s',           NumberValue(default=100, precision=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('samples_per_point',      NumberValue(default=50, precision=0, step=10, min=15, max=500), group='timing')

        # qubit
        self.setattr_argument('time_qubit_us',          NumberValue(default=5., precision=3, step=10, min=0.1, max=100000), group='qubit')
        self.setattr_argument("freq_qubit_mhz",         NumberValue(default=101.1038, precision=5, step=1, min=1, max=10000), group='qubit')
        self.setattr_argument("att_qubit_db",           NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5), group='qubit')

        # instantiate subsequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.rabiflop_subsequence =     RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.readout_subsequence =      ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        # get max readout time (b/c ReadoutAdaptive doesn't store it)
        time_readout_us =   self.get_parameter('time_readout_us', group='timing', override=False,
                                               conversion_function=us_to_mu)
        # calculate time taken for a single point
        self.time_per_point_us =    ((self.initialize_subsequence.time_repump_qubit_mu
                                     + self.initialize_subsequence.time_doppler_cooling_mu
                                     + self.initialize_subsequence.time_spinpol_mu
                                     + time_readout_us) / 1000
                                     + self.time_qubit_us)

        # get relevant timings and calculate the number of repetitions
        self.repetitions = round(self.time_total_s / (self.samples_per_point * self.time_per_point_us * us))

        # declare loop iterators and holder variables ahead of time to reduce overhead
        self._iter_repetitions =    np.arange(self.repetitions)
        self._iter_loop =           np.arange(self.samples_per_point)
        self._state_array =         np.zeros(self.samples_per_point, dtype=np.int32)

        # convert qubit parameters
        self.freq_qubit_ftw =   hz_to_ftw(self.freq_qubit_mhz * MHz)
        self.att_qubit_mu =     att_to_mu(self.att_qubit_db * dB)

    @rpc
    def initialize_plotting(self) -> TNone:
        """
        Configure datasets and applets for plotting.
        """
        # prepare datasets for storing counts
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        state_x_arr = np.zeros(self.repetitions) * np.nan
        state_x_arr[0] = 0
        state_y_arr = np.zeros(self.repetitions) * np.nan
        state_y_arr[0] = 0

        # prepare datasets for storing counts
        self.set_dataset('temp.qubit_align.counts_x', state_x_arr, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.qubit_align.counts_y', state_y_arr, broadcast=True, persist=False, archive=False)

        # initialize plotting applet
        self.ccb.issue(
            # name of broadcast & applet name
            "create_applet", "qubit_alignment",
            # command
            '${artiq_applet}plot_xy temp.qubit_align.counts_y'
            ' --x temp.qubit_align.counts_x --title "Qubit Alignment"',
            group=["alignment"] # folder directory for applet
        )

    @property
    def results_shape(self):
        return (self.repetitions, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # prepare plotting
        # note: do it here instead of prepare to prevent overriding other experiments
        self.initialize_plotting()
        self.core.break_realtime()

        # prepare qubit beam for readout
        self.qubit.set_profile(0)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit.set_mu(self.freq_qubit_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)
        delay_mu(10000)

        # record alignment sequence
        with self.core_dma.record('_QUBIT_ALIGNMENT'):
            # initialize ion in S-1/2 state & rabiflop
            self.initialize_subsequence.run()
            self.rabiflop_subsequence.run()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        # retrieve DMA handles for qubit alignment
        _handle_alignment = self.core_dma.get_handle('_QUBIT_ALIGNMENT')
        self.core.break_realtime()

        # predeclare variables
        ion_state = (-1, 0, np.int64(0))

        # MAIN LOOP
        for num_rep in self._iter_repetitions:
            delay_mu(25000)

            # average readings over N samples
            for num_count in self._iter_loop:
                delay_mu(25000)

                # run qubit alignment sequence
                self.core_dma.playback_handle(_handle_alignment)

                # determine ion state
                ion_state = self.readout_subsequence.run()
                self._state_array[num_count] = ion_state[0]
                delay_mu(25000)

            # update dataset
            self.update_results(num_rep, self._state_array)
            self.core.break_realtime()

            # periodically check termination
            if (num_rep % 10) == 0:
                self.check_termination()
                self.core.break_realtime()

    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num: TInt32, state_array: TArray(TInt32, 1)) -> TNone:
        """
        todo: document
        """
        # average results while ignoring indeterminate results
        num_indeterminate = np.sum(state_array[state_array == -1])
        if num_indeterminate == self.samples_per_point:
            dstate_probability = 0.
        else:
            dstate_probability = np.sum(state_array[state_array == 1]) / (self.samples_per_point - num_indeterminate)

        # update datasets for broadcast
        self.mutate_dataset('temp.qubit_align.counts_x', self._result_iter, iter_num * (self.samples_per_point * self.time_per_point_us * us))
        self.mutate_dataset('temp.qubit_align.counts_y', self._result_iter, dstate_probability)

        # update dataset for HDF5 storage
        self.mutate_dataset('results', self._result_iter,
                            np.array([iter_num * (self.samples_per_point * self.time_per_point_us * us),
                                      dstate_probability]))

        # update completion monitor
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._result_iter / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += 1

