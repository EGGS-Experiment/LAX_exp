import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon

# automatic thresholding
# qubit parameters: time, freq, power
# need to import doppler cool, rabiflop, readout


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
        self.setattr_argument('time_total_s',           NumberValue(default=20, ndecimals=0, step=100, min=5, max=100000), group='timing')
        self.setattr_argument('update_interval_ms',     NumberValue(default=500, ndecimals=0, step=100, min=50, max=100000), group='timing')
        self.setattr_argument('samples_per_point',      NumberValue(default=20, ndecimals=0, step=10, min=1, max=500), group='timing')

        # qubit
        self.setattr_argument('time_qubit_us',          NumberValue(default=3000, ndecimals=0, step=500, min=100, max=100000), group='qubit')
        self.setattr_argument("freq_qubit_mhz",         NumberValue(default=103.3455, ndecimals=5, step=1, min=1, max=10000), group='qubit')
        self.setattr_argument("att_qubit_db",           NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group='qubit')
        self.setattr_argument('counts_threshold',       NumberValue(default=60, ndecimals=0, step=10, min=1, max=1000), group='qubit')


        # instantiate subsequences
        self.initialize_subsequence =                   InitializeQubit(self)
        self.readout_subsequence =                      Readout(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        # ensure sample rates and readout times are reasonable
        # this is also necessary to ensure we have enough slack
        time_per_point_s =                      (2. * self.time_sample_us * us) * self.samples_per_point
        max_sampling_time_s =                   0.75 * (self.update_interval_ms * ms)
        assert time_per_point_s < max_sampling_time_s, "Error: unreasonable sample values idk - try again welp"
        # todo: also assert whether the time_slack_mu is valid

        # convert qubit parameters
        self.time_qubit_mu =                    self.core.seconds_to_mu(self.time_qubit_us * us)
        self.freq_qubit_ftw =                   hz_to_ftw(self.freq_qubit_mhz * MHz)
        self.att_qubit_mu =                     att_to_mu(self.att_qubit_db * dB)


        # get relevant timings and calculate the number of repetitions
        self.repetitions =                      round(self.time_total_s / (self.update_interval_ms * ms))
        # assume 50% duty cycle
        self.time_slack_mu =                    self.core.seconds_to_mu((self.update_interval_ms * ms) / self.samples_per_point
                                                                        - (2. * self.time_sample_us * us))

        # create holder variable that can store averaged counts
        self._counts_signal =                   np.int32(0)
        self._counts_background =               np.int32(0)

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
        return (self.repetitions, 2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # prepare qubit beam for readout
        self.qubit.set_profile(0)
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit.set_mu(self.freq_qubit_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)

        # record alignment sequence
        with self.core_dma.record('_QUBIT_ALIGNMENT'):
            # initialize ion in S-1/2 state
            self.initialize_subsequence.run()

            # run sideband cooling
            self.cooling_subsequence.run()

            # rabi flop
            self.qubit.on()
            delay_mu(self.time_qubit_mu)
            self.qubit.off()

            # readout
            self.readout_subsequence.run()
        
        
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()

        # retrieve DMA handles for qubit alignment
        _handle_alignment = self.core_dma.get_handle('_QUBIT_ALIGNMENT')
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


    # overload the update_results function to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num, counts):
        self.mutate_dataset('_tmp_counts_x', self._result_iter, iter_num * (self.update_interval_ms * ms))
        self.mutate_dataset('_tmp_counts_y', self._result_iter, np.array([counts / self.samples_per_point]))
        self.set_dataset('management.completion_pct', round(100. * self._result_iter / len(self.results), 3), broadcast=True, persist=True, archive=False)
        self._result_iter += 1

    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
