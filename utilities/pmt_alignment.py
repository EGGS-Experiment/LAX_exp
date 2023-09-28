import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment


class PMTAlignment(LAXExperiment, Experiment):
    """
    Utility: PMT Alignment

    Read PMT counts over time with cooling repump on/off to compare signal/background.
    """
    name = 'PMT Alignment'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # timing
        self.setattr_argument('time_total_s',           NumberValue(default=1000, ndecimals=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('update_interval_ms',     NumberValue(default=500, ndecimals=0, step=100, min=50, max=100000), group='timing')

        self.setattr_argument('time_sample_us',         NumberValue(default=3000, ndecimals=0, step=500, min=100, max=100000), group='timing')
        self.setattr_argument('samples_per_point',      NumberValue(default=50, ndecimals=0, step=10, min=1, max=500), group='timing')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        # ensure sample rates and readout times are reasonable
        # this is also necessary to ensure we have enough slack
        time_per_point_s =                          (2. * self.time_sample_us * us) * self.samples_per_point
        max_sampling_time_s =                       0.75 * (self.update_interval_ms * ms)
        assert time_per_point_s < max_sampling_time_s, "Error: unreasonable sample values idk - try again welp"
        # todo: also assert whether the time_slack_mu is valid

        # get relevant timings and calculate the number of repetitions
        self.repetitions =                      round(self.time_total_s / (self.update_interval_ms * ms))
        self.time_sample_mu =                   self.core.seconds_to_mu(self.time_sample_us * us)
        # assume 50% duty cycle
        self.time_slack_mu =                    self.core.seconds_to_mu((self.update_interval_ms * ms) / self.samples_per_point
                                                                        - (2. * self.time_sample_us * us))

        # create holder variable that can store averaged counts
        self._counts_signal =                   np.int32(0)
        self._counts_background =               np.int32(0)

        # tmp remove
        # declare the loop iterators ahead of time to reduce overhead
        self._iter_repetitions =                np.arange(self.repetitions)
        self._iter_loop =                       np.arange(self.samples_per_point)

        # prepare datasets for storing counts
        self.set_dataset('_tmp_counts_x',     np.zeros(self.repetitions), broadcast=True, persist=False, archive=False)
        self.set_dataset('_tmp_counts_y',     np.zeros((self.repetitions, 3)), broadcast=True, persist=False, archive=False)
        self.setattr_dataset('_tmp_counts_x')
        self.setattr_dataset('_tmp_counts_y')


    @property
    def results_shape(self):
        return (self.repetitions, 4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # ensure beams are set to readout
        # todo: more rigorous setting/checking of beam waveforms etc.
        self.pump.readout()

        # record alignment sequence
        with self.core_dma.record('_PMT_ALIGNMENT'):
            # todo: should we set ensure readout beams are set here?

            # set beams to record signal + background
            self.repump_cooling.on()
            self.pump.on()

            # record PMT signal + background
            self.pmt.count(self.time_sample_mu)

            # set beams to record signal + background
            self.repump_cooling.off()

            # add holdoff time to allow ions/beams/idk to settle
            delay_mu(2500)

            # record PMT background only
            self.pmt.count(self.time_sample_mu)

            # turn repump beams back on to cool ions
            self.repump_cooling.on()
        
        
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()

        # retrieve DMA handles for PMT alignment
        _handle_alignment = self.core_dma.get_handle('_PMT_ALIGNMENT')
        self.core.break_realtime()


        # MAIN LOOP
        for i in self._iter_repetitions:

            # clear holder variables
            self._counts_signal =       0
            self._counts_background =   0
            self.core.break_realtime()

            # average readings over N samples
            for j in self._iter_loop:

                # run PMT alignment sequence
                self.core_dma.playback_handle(_handle_alignment)

                # fetch running totals of PMT counts
                self._counts_signal +=      self.pmt.fetch_count()
                self._counts_background +=  self.pmt.fetch_count()

                # wait until next sample
                delay_mu(self.time_slack_mu)

            # update dataset
            with parallel:
                self.update_results(i, self._counts_signal, self._counts_background)
                self.core.break_realtime()

    # tmp remove
    @rpc(flags={"async"})
    def update_results(self, iter_num, counts_signal, counts_background):
        self.mutate_dataset('_tmp_counts_x', self._result_iter, iter_num * (self.update_interval_ms * ms))
        self.mutate_dataset('_tmp_counts_y', self._result_iter, np.array([counts_signal / self.samples_per_point,
                                                                          counts_background / self.samples_per_point,
                                                                          (counts_signal - counts_background) / self.samples_per_point]))
        self.set_dataset('management.completion_pct', round(100. * self._result_iter / len(self.results), 3), broadcast=True, persist=True, archive=False)
        self._result_iter += 1

    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # print results
        # print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
