import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment


# todo: sequence
# *** check wavemeter freqs all OK
# *** set low voltages
# *** set cooling beams & turn off b-field (to use spinpol)
# *** start timer & oven
# *** read PMT counts until they start to go high for a while
# *** stop oven/increase electrodes/stop timer
# *** reduce beam powers
# ** integrate labjack


class LoadIon(LAXExperiment, Experiment):
    """
    Utility: Load Ion

    Automated loading of an ion via the Calcium oven.
    """
    name = 'Load Ion'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # loading
        self.setattr_argument('time_total_s',           NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='loading')
        self.setattr_argument('monitor_interval_ms',    NumberValue(default=500, ndecimals=0, step=100, min=50, max=100000), group='loading')
        self.setattr_argument('oven_current_amps',      NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='loading')
        self.setattr_argument('ion_threshold_counts',   NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='loading')

        # electrode voltages
        self.setattr_argument('voltage_load_east_endcap_v', NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_loading')
        self.setattr_argument('voltage_load_west_endcap_v', NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_loading')
        self.setattr_argument('voltage_load_vert_shim_v',   NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_loading')
        self.setattr_argument('voltage_load_horiz_shim_v',  NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_loading')
        self.setattr_argument('voltage_load_a_ramp_v',      NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_loading')

        # trapping voltages
        self.setattr_argument('voltage_trap_east_endcap_v', NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_trapping')
        self.setattr_argument('voltage_trap_west_endcap_v', NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_trapping')
        self.setattr_argument('voltage_trap_vert_shim_v',   NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_trapping')
        self.setattr_argument('voltage_trap_horiz_shim_v',  NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_trapping')
        self.setattr_argument('voltage_trap_a_ramp_v',      NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='voltages_trapping')

        # relevant devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')


    def prepare_experiment(self):
        # convert relevant timings to machine units
        self.time_slack_mu =        self.core.seconds_to_mu(24 * us)
        self.time_sample_mu =       self.core.seconds_to_mu(self.time_sample_ms * ms)
        self.time_per_point_mu =    (self.time_sample_ms + self.time_slack_mu) * (self.signal_samples_per_point + self.bgr_samples_per_point)

        # calculate the number of repetitions
        self.repetitions =          round(self.time_total_s / self.core.mu_to_seconds(self.time_per_point_mu))

        # create holder variable that can store averaged counts
        self._counts_signal =       np.int32(0)
        self._counts_background =   np.int32(0)

        # declare the loop iterators ahead of time to reduce overhead
        self._iter_repetitions =    np.arange(self.repetitions)
        self._iter_signal =         np.arange(self.signal_samples_per_point)
        self._iter_background =     np.arange(self.background_samples_per_point)

        # prepare datasets for storing counts
        self.set_dataset('temp.imag_align._tmp_counts_x',  np.zeros(self.repetitions), broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.imag_align._tmp_counts_y',  np.zeros((self.repetitions, 3)), broadcast=True, persist=False, archive=False)

    @property
    def results_shape(self):
        return (self.repetitions, 4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # configure beams
        self.pump.readout()
        self.pump.on()
        self.repump_qubit.on()

        # record alignment sequence - signal
        with self.core_dma.record('_PMT_ALIGNMENT_SIGNAL'):
            # activate doppler repump beam
            self.repump_cooling.on()

            # get signal counts
            for i in self._iter_signal:
                self.pmt.count(self.time_sample_mu)
                self.pump.off()


        # record alignment sequence - background
        with self.core_dma.record('_PMT_ALIGNMENT_BACKGROUND'):
            # disable doppler repump beam
            self.repump_cooling.off()

            # get background counts
            for i in self._iter_background:
                self.pmt.count(self.time_sample_mu)
                self.pump.off()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()

        # retrieve DMA handles for PMT alignment
        _handle_alignment_signal = self.core_dma.get_handle('_PMT_ALIGNMENT_SIGNAL')
        _handle_alignment_background = self.core_dma.get_handle('_PMT_ALIGNMENT_BACKGROUND')
        self.core.break_realtime()


        # MAIN LOOP
        for i in self._iter_repetitions:

            # clear holder variables
            self._counts_signal =       0
            self._counts_background =   0
            self.core.break_realtime()

            # store signal counts
            self.core_dma.playback_handle(_handle_alignment_signal)
            # retrieve signal counts
            for i in self._iter_signal:
                self._counts_signal += self.pmt.fetch_count()

            # store background counts
            self.core_dma.playback_handle(_handle_alignment_background)
            # retrieve background counts
            for i in self._iter_background:
                self._counts_background += self.pmt.fetch_count()

            # # add slack
            # delay_mu(self.time_slack_mu)

            # update dataset
            with parallel:
                self.update_results(i, self._counts_signal, self._counts_background)
                self.core.break_realtime()


    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num, counts_signal, counts_background):
        self.mutate_dataset('_tmp_counts_x', self._result_iter, iter_num * (self.update_interval_ms * ms))
        self.mutate_dataset('_tmp_counts_y', self._result_iter, np.array([counts_signal / self.samples_per_point,
                                                                          counts_background / self.samples_per_point,
                                                                          (counts_signal - counts_background) / self.samples_per_point]))
        self.set_dataset('management.dynamic.completion_pct', round(100. * self._result_iter / len(self.results), 3), broadcast=True, persist=True, archive=False)
        self._result_iter += 1



    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # print results
        # print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
