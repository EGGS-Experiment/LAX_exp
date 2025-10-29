import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt
# todo: create a slower, vibes-based sampler read exp
# todo: implement progress bar


class SamplerReadFast(EnvExperiment):
    """
    Tool: Sampler Read Fast

    Read Sampler values over time.
    Warning: due to timing concerns, this experiment CAN NOT be cancelled after submission
        without loss of data.
    """
    kernel_invariants = {
        'adc', 'time_delay_mu', 'repetitions',
        'channel_list', 'channel_iter', 'gain_list_mu', 'adc_mu_to_v_list',
    }

    def build(self):
        # devices
        self.setattr_device('core')

        # channels & gains
        self.setattr_argument("channel_gain_dict",  PYONValue({2: 100}),
                              tooltip="A dictionary of {channel_number: gain_factor}.\n"
                                      "Gain factor must be an int in [1, 10, 100, 1000].\n"
                                      "Can record multiple channels, though this will limit the max sample rate (due to timing overheads).")

        # timing
        self.setattr_argument("sample_rate_hz", NumberValue(default=5000, precision=3, step=1, min=1, max=5100, scale=1., unit="Hz"),
                              tooltip="The ADC sample rate.\n"
                                      "This is limited to 5kHz for long-term measurement.\n"
                                      "For shorter timescales (e.g. ~1s), a sample rate of up to 20kHz can be achieved.")
        self.setattr_argument("time_total_s",   NumberValue(default=1, precision=0, step=1, min=1, max=100000, scale=1., unit="s"),
                              tooltip="The total length of time to record for.\n"
                                      "Note: once this experiment begins, it can not be cancelled without loss of data, "
                                      "and termination CANNOT be requested.\n"
                                      "If possible, a long run should be broken up into multiple smaller runs.")

    def prepare(self):
        # general
        self.channel_list =     list(self.channel_gain_dict.keys())
        self.channel_iter =     list(range(len(self.channel_gain_dict)))
        self.gain_list_mu =     [int(np.log10(gain_mu)) for gain_mu in self.channel_gain_dict.values()]
        self.adc_mu_to_v_list = np.array([10 / (2**15 * gain_mu) for gain_mu in self.channel_gain_dict.values()])

        # ADC
        self.adc = self.get_device("sampler0")

        # timing
        self.time_delay_mu =    self.core.seconds_to_mu(1 / self.sample_rate_hz)
        self.repetitions =      np.int32(self.time_total_s * self.sample_rate_hz)

        # datasets
        self.set_dataset('results', np.zeros([self.repetitions, len(self.channel_list)]))
        self.setattr_dataset('results')

        # save parameters
        self.set_dataset('sample_rate_hz', self.sample_rate_hz)
        self.set_dataset('time_total_s', self.time_total_s)

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # set ADC channel gains
        for i in self.channel_iter:
            self.adc.set_gain_mu(self.channel_list[i], self.gain_list_mu[i])
            self.core.break_realtime()

        # create holding buffer
        sampler_buffer = [0] * 8
        delay_mu(10000)

        # sampling loop
        for i in range(self.repetitions):
            with parallel:
                delay_mu(self.time_delay_mu)
                with sequential:
                    self.adc.sample_mu(sampler_buffer)
                    self.update_dataset(i, sampler_buffer)

        # ensure events finish submission before experiment completion
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())

    @rpc(flags={"async"})
    def update_dataset(self, i, volts_mu_arr):
        """
        Records values via rpc to minimize kernel overhead.
        :param i: todo: document
        :param volts_mu_arr: todo: document
        """
        data = np.array(volts_mu_arr)[self.channel_list] * self.adc_mu_to_v_list
        self.mutate_dataset("results", i, data)

    def analyze(self):
        # print out statistics of results
        print('\tResults:')
        for i in self.channel_iter:
            print('\t\tCH{:d}:\t{:.3f} +/- {:.3f} mV'.format(self.channel_list[i],
                                                             np.mean(self.results[:, i]) * 1000,
                                                             np.std(self.results[:, i]) * 1000))

