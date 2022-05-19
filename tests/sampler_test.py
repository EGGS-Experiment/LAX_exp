import numpy as np
from artiq.experiment import *


class sampler_exp(EnvExperiment):
    """
    Read out data from the Sampler.
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device('sampler0')

    def prepare(self):
        self.samples = 1000
        self.set_data("sampler_readout", np.zeros(self.samples), broadcast=True)
        self.readout_channel = 0

    @kernel
    def run(self):
        """
        todo
        """
        self.core.reset()
        # set up ADC
        self.sampler0.init()
        self.sampler0.set_gain_mu(self.readout_channel, 1)
        self.core.break_realtime()
        # create buffer
        sampler_buffer = [0] * 8
        for i in range(self.samples):
            self.sampler0.sample_mu(sampler_buffer)
            self.mutate_dataset("sampler_readout", i, sampler_buffer[self.readout_channel])
            self.core.break_realtime()

    def analyze(self):
        """
        todo:
        """
        # todo: convert from machine units to volts
        pass
