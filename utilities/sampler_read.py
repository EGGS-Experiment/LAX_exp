import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class sampler_exp(EnvExperiment):
    """
    Sampler Read
    Read out data from the Sampler.
    """

    def build(self):
        # devices
        self.setattr_device('core')
        self.setattr_device('sampler0')

        # arguments
        self.setattr_argument("channel_readout", NumberValue(default=1, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("channel_gain_10dB", NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("time_delay_us", NumberValue(default=500, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("num_samples", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=20000))


    def prepare(self):
        # ADC
        self.adc = self.sampler0
        self.adc_mu_to_volts = (10 ** (1 - self.channel_gain_10dB)) / (2 ** 15)

        # timing
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # datasets
        self.set_dataset('sampler_readout', np.zeros(self.num_samples), broadcast=True)
        self.setattr_dataset('sampler_readout')

    @kernel
    def run(self):
        self.core.reset()

        # set up ADC
        self.sampler0.init()
        self.sampler0.set_gain_mu(self.channel_readout, self.channel_gain_10dB)
        sampler_buffer = [0] * 8
        self.core.break_realtime()

        # sampling loop
        for i in range(self.num_samples):
            with parallel:
                delay_mu(self.time_delay_mu)
                self.sampler0.sample_mu(sampler_buffer)
            self.mutate_dataset("sampler_readout", i, sampler_buffer[self.channel_readout] * self.adc_mu_to_volts)

    def analyze(self):
        print("sampler ch{:d}: {:.4f} +/- {:.5f}".format(self.channel_readout, np.average(self.sampler_readout), np.std(self.sampler_readout)))
