import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class sampler_exp(EnvExperiment):
    """
    Read out data from the Sampler.
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device('sampler0')

    def prepare(self):
        # values
        self.samples = 50
        self.set_dataset("sampler_readout", np.zeros(self.samples), broadcast=True)
        self.setattr_dataset('sampler_readout')
        # sampler
        self.readout_channel = 0
        self.readout_gain = 1
        # readout interval
        self.time_delay_mu = self.core.seconds_to_mu(100 * us)

    @kernel
    def run(self):
        """
        todo
        """
        self.core.reset()
        # set up ADC
        self.sampler0.init()
        self.sampler0.set_gain_mu(self.readout_channel, self.readout_gain)
        self.core.break_realtime()
        # create buffer
        sampler_buffer = [0] * 8
        for i in range(self.samples):
            with parallel:
                delay_mu(self.time_delay_mu)
                with sequential:
                    self.sampler0.sample_mu(sampler_buffer)
                    self.core.break_realtime()
                    self.mutate_dataset("sampler_readout", i, sampler_buffer[self.readout_channel])

    def analyze(self):
        """
        Convert values from machine units to volts.
        """
        print('sampler readout:', self.sampler_readout)
        print('sampler avg:', np.average(self.sampler_readout)/2**15)
        # dataset = self.get_dataset("sampler_readout")
        # for i, val in enumerate(dataset):
        #     self.mutate_dataset("sampler_readout", i, adc_mu_to_volt(val))
