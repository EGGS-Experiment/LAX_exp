from artiq.experiment import *
from numpy import zeros, int32, arange, array, nan
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment


class ImagingAlignment(LAXExperiment, Experiment):
    """
    Utility: Imaging Alignment

    Read PMT counts over time with cooling repump on/off to compare signal/background.
    """
    name = 'Imaging Alignment'
    kernel_invariants = {
        # conversions, counters etc.
        "repetitions", "_iter_repetitions", "_iter_signal", "_iter_background",

        # timing
        "time_slack_us", "time_per_point_s", "time_per_point_mu", "time_slack_mu", "time_sample_mu",

        # other hardware values
        "freq_readout_ftw", "ampl_readout_asf"
    }

    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # general
        self.setattr_argument('time_total_s', NumberValue(default=10, precision=0, step=100, min=5, max=100000, scale=1., unit="s"))

        # sampling
        self.setattr_argument('signal_samples_per_point',       NumberValue(default=48, precision=0, step=10, min=1, max=100),
                              group='sampling')
        self.setattr_argument('background_samples_per_point',   NumberValue(default=5, precision=0, step=2, min=1, max=100),
                              group='sampling')

        # readout
        self.setattr_argument('time_sample_us',     NumberValue(default=3000, precision=1, step=500, min=100, max=100000, scale=1., unit="us"),
                              group='readout')
        self.setattr_argument("freq_readout_mhz",   NumberValue(default=108., precision=6, step=1, min=1, max=500, scale=1., unit='MHz'),
                              group='cooling')
        self.setattr_argument("ampl_readout_pct",   NumberValue(default=42., precision=2, step=5, min=0.01, max=50, scale=1., unit='%'),
                              group='cooling')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        # convert relevant values to machine units
        self.time_slack_us =    2   # determined empirically
        self.time_per_point_s = ((self.time_slack_us * us + self.time_sample_us * us) *
                                 (self.signal_samples_per_point + self.background_samples_per_point))

        self.time_slack_mu =        self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_sample_mu =       self.core.seconds_to_mu(self.time_sample_us * us)
        self.time_per_point_mu =    self.core.seconds_to_mu(self.time_per_point_s)

        self.freq_readout_ftw = self.pump.frequency_to_ftw(self.freq_readout_mhz * MHz)
        self.ampl_readout_asf = self.pump.frequency_to_ftw(self.ampl_readout_pct / 100.)

        # calculate number of repetitions
        self.repetitions =          round(self.time_total_s / self.core.mu_to_seconds(self.time_per_point_mu))

        # create holder variable to store averaged counts
        self._counts_signal =       int32(0)
        self._counts_background =   int32(0)

        # declare the loop iterators ahead of time to reduce overhead
        self._iter_repetitions =    arange(self.repetitions)
        self._iter_signal =         arange(self.signal_samples_per_point)
        self._iter_background =     arange(self.background_samples_per_point)

    @rpc
    def initialize_plotting(self) -> TNone:
        """
        Configure datasets and applets for plotting.
        """
        # prepare datasets for storing counts
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        counts_x_arr = zeros(self.repetitions) * nan
        counts_x_arr[0] = 0
        counts_y_arr = zeros((self.repetitions, 3)) * nan
        counts_y_arr[0, :] = 0
        counts_snr_arr = zeros(self.repetitions) * nan
        counts_snr_arr[0] = 0

        self.set_dataset('temp.imag_align.counts_x', counts_x_arr, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.imag_align.counts_y', counts_y_arr, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.imag_align.counts_snr', counts_snr_arr, broadcast=True, persist=False, archive=False)

        # initialize plotting applet
        self.ccb.issue(
            "create_applet",        # name of broadcast
            "imaging_alignment",    # applet name
            # command
            '$python -m LAX_exp.applets.plot_xy_multi temp.imag_align.counts_y'
            ' --x temp.imag_align.counts_x --title "Imaging Alignment"',
            group=["alignment"] # folder directory for applet
        )
        self.ccb.issue(
            "create_applet",            # name of broadcast
            "imaging_alignment_SNR",    # applet name
            # command
            '${artiq_applet}plot_xy temp.imag_align.counts_snr'
            ' --x temp.imag_align.counts_x --title "Imaging Alignment - SNR"',
            group=["alignment"] # folder directory for applet
        )

    @property
    def results_shape(self):
        return (self.repetitions, 4)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # prepare plotting
        # note: do it here instead of prepare to prevent overriding other experiments
        self.initialize_plotting()
        self.core.break_realtime()

        # configure beams
        self.pump.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf,
                         profile=self.pump.profile_readout,
                         phase_mode=PHASE_MODE_CONTINUOUS)
        self.pump.readout()
        self.pump.on()
        self.repump_qubit.on()

        # record alignment sequence - signal
        with self.core_dma.record('_PMT_ALIGNMENT_SIGNAL'):
            # activate doppler repump beam
            self.repump_cooling.on()

            # get signal counts
            for _ in self._iter_signal:
                self.pmt.count(self.time_sample_mu)
                delay_mu(self.time_slack_mu)

        # record alignment sequence - background
        with self.core_dma.record('_PMT_ALIGNMENT_BACKGROUND'):
            # disable doppler repump beam
            self.repump_cooling.off()

            # get background counts
            for _ in self._iter_background:
                self.pmt.count(self.time_sample_mu)
                delay_mu(self.time_slack_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # retrieve DMA handles for PMT alignment sequences
        _handle_alignment_signal =      self.core_dma.get_handle('_PMT_ALIGNMENT_SIGNAL')
        _handle_alignment_background =  self.core_dma.get_handle('_PMT_ALIGNMENT_BACKGROUND')

        # MAIN LOOP
        for num_rep in self._iter_repetitions:
            # clear holder variables
            self._counts_signal =       0
            self._counts_background =   0
            self.core.break_realtime()

            # get signal
            self.core_dma.playback_handle(_handle_alignment_signal)
            # retrieve signal counts
            for _ in self._iter_signal:
                self._counts_signal += self.pmt.fetch_count()
            self.core.break_realtime()

            # get background
            self.core_dma.playback_handle(_handle_alignment_background)
            # retrieve background counts
            for _ in self._iter_background:
                self._counts_background += self.pmt.fetch_count()

            # update dataset & periodically check termination
            self.update_results(num_rep, self._counts_signal, self._counts_background)
            if (num_rep % 10) == 0:
                self.check_termination()

    # overload update_results to allow real-time dataset updating
    @rpc(flags={"async"})
    def update_results(self, iter_num: TInt32, counts_signal: TInt64, counts_background: TInt64) -> TNone:
        # convert total counts into averaged counts
        _counts_avg_signal =        counts_signal / self.signal_samples_per_point
        _counts_avg_background =    counts_background / self.background_samples_per_point

        # update datasets for broadcast
        self.mutate_dataset('temp.imag_align.counts_x', self._result_iter, iter_num * self.time_per_point_s)
        self.mutate_dataset('temp.imag_align.counts_y', self._result_iter, array([_counts_avg_signal,
                                                                                  _counts_avg_background,
                                                                                  _counts_avg_signal - _counts_avg_background]))
        self.mutate_dataset('temp.imag_align.counts_snr',
                            self._result_iter, (_counts_avg_signal - _counts_avg_background) / _counts_avg_background)

        # update dataset for HDF5 storage
        self.mutate_dataset('results', self._result_iter,
                            array([iter_num * self.time_per_point_s,
                                      _counts_avg_signal,
                                      _counts_avg_background,
                                      _counts_avg_signal - _counts_avg_background]))

        # update completion monitor
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._result_iter / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._result_iter += 1

