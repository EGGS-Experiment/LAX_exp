import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment


class ImagingAlignment(LAXExperiment, Experiment):
    """
    Utility: Imaging Alignment

    Read PMT counts over time with cooling repump on/off to compare signal/background.
    """
    name = 'Imaging Alignment'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # timing
        self.setattr_argument('time_total_s',               NumberValue(default=200, precision=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('time_sample_us',             NumberValue(default=3000, precision=1, step=500, min=100, max=100000), group='timing')

        # sampling
        self.setattr_argument('signal_samples_per_point',       NumberValue(default=48, precision=0, step=10, min=1, max=100), group='sampling')
        self.setattr_argument('background_samples_per_point',   NumberValue(default=5, precision=0, step=2, min=1, max=100), group='sampling')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        # convert relevant timings to machine units
        self.time_slack_us =        2
        self.time_per_point_s =     (self.time_slack_us * us + self.time_sample_us * us)* (self.signal_samples_per_point + self.background_samples_per_point)

        self.time_slack_mu =        self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_sample_mu =       self.core.seconds_to_mu(self.time_sample_us * us)
        self.time_per_point_mu =    self.core.seconds_to_mu(self.time_per_point_s)

        # calculate the number of repetitions
        self.repetitions =          round(self.time_total_s / self.core.mu_to_seconds(self.time_per_point_mu))

        # create holder variable that can store averaged counts
        self._counts_signal =       np.int32(0)
        self._counts_background =   np.int32(0)

        # declare the loop iterators ahead of time to reduce overhead
        self._iter_repetitions =    np.arange(self.repetitions)
        self._iter_signal =         np.arange(self.signal_samples_per_point)
        self._iter_background =     np.arange(self.background_samples_per_point)

        # todo: move to prepare/run somehow to prevent overriding
        # prepare datasets for storing counts
        self.set_dataset('temp.imag_align.counts_x', np.zeros(self.repetitions) * np.nan, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.imag_align.counts_y', np.zeros((self.repetitions, 3)) * np.nan, broadcast=True, persist=False, archive=False)


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
                delay_mu(self.time_slack_mu)


        # record alignment sequence - background
        with self.core_dma.record('_PMT_ALIGNMENT_BACKGROUND'):
            # disable doppler repump beam
            self.repump_cooling.off()

            # get background counts
            for i in self._iter_background:
                self.pmt.count(self.time_sample_mu)
                delay_mu(self.time_slack_mu)

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
            for j in self._iter_signal:
                self._counts_signal += self.pmt.fetch_count()
            self.core.break_realtime()

            # store background counts
            self.core_dma.playback_handle(_handle_alignment_background)
            # retrieve background counts
            for j in self._iter_background:
                self._counts_background += self.pmt.fetch_count()
            self.core.break_realtime()

            # # add slack
            # delay_mu(self.time_slack_mu)

            # update dataset
            with parallel:
                self.update_results(i, self._counts_signal, self._counts_background)
                self.core.break_realtime()


    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num: TInt32, counts_signal: TInt64, counts_background: TInt64) -> TNone:
        # convert total counts into averaged counts
        _counts_avg_signal =        counts_signal / self.signal_samples_per_point
        _counts_avg_background =    counts_background / self.background_samples_per_point

        # update datasets
        self.mutate_dataset('temp.imag_align.counts_x', self._result_iter, iter_num * self.time_per_point_s)
        self.mutate_dataset('temp.imag_align.counts_y', self._result_iter, np.array([_counts_avg_signal,
                                                                                          _counts_avg_background,
                                                                                          _counts_avg_signal - _counts_avg_background]))

        # update completion monitor
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._result_iter / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += 1


    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        print(self.get_dataset('temp.imag_align._tmp_counts_x'))
        # print results
        # print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
