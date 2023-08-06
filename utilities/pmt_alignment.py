from artiq.experiment import *
from numpy import zeros, arange, mean, std, int32, array

_DMA_HANDLE = 'PMT_ALIGN'


class pmt_alignment(EnvExperiment):
    """
    PMT Alignment

    Read PMT counts over time with cooling repump on/off to compare signal/background.
    """
    kernel_invariants = {
        'time_bin_mu',
        'time_reset_mu'
    }

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        # set necessary devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("urukul1_cpld")

        # timing
        self.setattr_argument('time_total_s',                       NumberValue(default=10000, ndecimals=6, step=1, min=0, max=100000), group='timing')
        self.setattr_argument('time_bin_us',                        NumberValue(default=500, ndecimals=3, step=1, min=0.01, max=10000000), group='timing')
        self.setattr_argument("sample_rate_hz",                     NumberValue(default=1000, ndecimals=3, step=1, min=1, max=100000), group='timing')
        # todo: create update iter argument which sets update rate

        # PMT
        self.setattr_argument("pmt_input_channel",                  NumberValue(default=0, ndecimals=0, step=1, min=0, max=7), group='pmt')
        self.setattr_argument("pmt_gating_edge",                    EnumerationValue(["rising", "falling", "both"], default="rising"), group='pmt')


    def prepare(self):
        """
        Set up the dataset and prepare things such that the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # iterators
        self.loop_iter =                                            arange(self.time_total_s * self.sample_rate_hz, dtype=int32)
        self.update_iter =                                          arange(self.time_total_s * self.sample_rate_hz, dtype=int32)

        # timing
        self.time_bin_mu =                                          self.core.seconds_to_mu(self.time_bin_us * us)
        self.time_reset_mu =                                        self.core.seconds_to_mu(1 / self.sample_rate_hz - self.time_bin_us * us)

        # set up datasets
        self.set_dataset('pmt_dataset',                             zeros((len(self.loop_iter)), 3))
        self.setattr_dataset('pmt_dataset')

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # program pulse sequence onto core DMA
        self.DMArecord()
        self.core.break_realtime()

        # retrieve pulse sequence handle
        handle = self.core_dma.get_handle(_DMA_HANDLE)
        self.core.break_realtime()

        # run the experiment
        for i in self.loop_iter:

            # run pulse sequence from core DMA
            self.core_dma.playback_handle(handle)

            # record pmt counts to dataset
            with parallel:
                self.update_dataset(i, self.pmt_counter.fetch_count(), self.pmt_counter.fetch_count())
                delay_mu(self.time_reset_mu)

    @kernel
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE):
            # record total signal
            self.urukul1_cpld.cfg_switches(0b1110)
            self.pmt_gating_edge(self.time_bin_mu)

            # turn off repump to record background signal
            self.urukul1_cpld.cfg_switches(0b0110)
            self.pmt_gating_edge(self.time_bin_mu)

    @rpc(flags={"async"})
    def updateDataset(self, index, counts_signal, counts_background):
        # create update array todo: document better
        update_arr = array([counts_signal, counts_background, counts_signal - counts_background])

        # update dataset
        self.mutate_dataset('pmt_dataset', index, update_arr)

        # todo: if index % update_iter, update labrad


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # print results
        print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
