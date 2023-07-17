from artiq.experiment import *
from numpy import zeros, arange, mean, std, int32

_DMA_HANDLE = 'PMT_exp'


class counter_read(EnvExperiment):
    """
    Counter Read

    Read TTL counts over time.
    """
    kernel_invariants = {
        'time_bin_mu',
        'time_reset_mu'
    }

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # timing
        self.setattr_argument('time_total_s',                       NumberValue(default=10000, ndecimals=6, step=1, min=0, max=100000))
        self.setattr_argument('time_bin_us',                        NumberValue(default=500, ndecimals=3, step=1, min=0.01, max=10000000))
        self.setattr_argument("sample_rate_hz",                     NumberValue(default=1000, ndecimals=3, step=1, min=1, max=100000))

        # PMT
        self.setattr_argument("pmt_input_channel",                  NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("pmt_gating_edge",                    EnumerationValue(["rising", "falling", "both"], default="rising"))


    def prepare(self):
        """
        Set up the dataset and prepare things such that the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # iterators
        self.loop_iter =                                            arange(self.time_total_s * self.sample_rate_hz, dtype=int32)

        # timing
        self.time_bin_mu =                                          self.core.seconds_to_mu(self.time_bin_us * us)
        self.time_reset_mu =                                        self.core.seconds_to_mu(1 / self.sample_rate_hz - self.time_bin_us * us)

        # set up datasets
        self.set_dataset('pmt_dataset',                             zeros(len(self.loop_iter)))
        self.setattr_dataset('pmt_dataset')

    @kernel
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
                self.mutate_dataset('pmt_dataset', i, self.pmt_counter.fetch_count())
                delay_mu(self.time_reset_mu)

    @kernel
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE):
            self.pmt_gating_edge(self.time_bin_mu)


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # print results
        print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
