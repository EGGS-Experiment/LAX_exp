import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class sampler_fast(EnvExperiment):
    """
    Sampler Fast
    Read out data from the Sampler as fast as we can.
    """

    def build(self):
        # devices
        self.setattr_device('core')
        self.setattr_device('sampler0')

        # arguments
        self.setattr_argument("channel_readout", NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("channel_gain_10dB", NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("time_delay_us", NumberValue(default=7, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("num_samples", NumberValue(default=10000, ndecimals=0, step=1, min=1, max=20000))

    def prepare(self):
        # ADC
        self.adc = self.sampler0
        self.adc_buffer = [0] * 2
        self.adc_mu_to_volts = (10 ** (1 - self.channel_gain_10dB)) / (2 ** 15)

        # timing
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # datasets
        self.set_dataset('sampler_readout', [], broadcast=True)
        self.setattr_dataset('sampler_readout')

        # iterables
        self.range = list(range(self.num_samples))

    @kernel(flags="fast-math")
    def run(self):
        self.core.reset()

        # set up ADC
        self.adc.init()
        self.adc.set_gain_mu(self.channel_readout, self.channel_gain_10dB)
        self.core.break_realtime()

        # sampling loop
        # for _ in self.range:
        #for _ in self.range:
        while True:
            with parallel:
                delay_mu(self.time_delay_mu)
                #with sequential:
                self.sampler0.sample_mu_fast(self.adc_buffer)
                    #self.update_dataset(self.adc_buffer[self.channel_readout])

    @rpc(flags={"async"})
    def update_dataset(self, adc_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset("sampler_readout", self.adc_mu_to_volts)

    def analyze(self):
        pass
        #print("sampler ch{:d}: {:f} +/- {:f}".format(self.channel_readout, np.average(self.sampler_readout), np.std(self.sampler_readout)))
