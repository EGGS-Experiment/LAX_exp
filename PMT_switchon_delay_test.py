from artiq.experiment import *
from numpy import int32, int64
import numpy as np


class PMTSwitchonTest(EnvExperiment):
    """
    PMT Switch On Test.
    """

    def build(self):
        #self.setattr_argument("frequency", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("amplitude", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("phase", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("channel", NumberValue(ndecimals=0, step=1))
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl1")
        self.setattr_device("ttl4")
        self.setattr_device("ttl5")
        self.num_bins = 1000
        self.time_step = self.core.seconds_to_mu(1000 * ns)
        self.switchon_delay = self.core.seconds_to_mu(1000 * ns)
        #self.reset_time = 1 * ms
        #ttl0 count read
        #ttl1 over light
        #ttl2 linetrigger in
        #ttl4 power
        #ttl5 linetrigger out

    def prepare(self):
        #self.time_step_mu = self.core.seconds_to_mu(1 * us)
        # set dataset
        self.dataset_name = 'oktc'
        #self.set_dataset(self.dataset_name, np.zeros(self.num_bins), broadcast=True)
        self.set_dataset(self.dataset_name, np.zeros(self.num_bins), broadcast=True)
        self.setattr_dataset(self.dataset_name)

    @kernel
    def run(self):
        self.core.reset()
        #self.ttl0.input()
        #self.core.reset()
        for i in range(self.num_bins):
            self.core.break_realtime()
            self.ttl4.on()
            #self.core.reset()
            delay_mu(i * self.switchon_delay)
            self.mutate_dataset(self.dataset_name, i, self.ttl0.count(self.ttl0.gate_rising_mu(self.time_step)))
            self.core.break_realtime()
            self.ttl4.off()
            self.core.break_realtime()
            # get over light
            #self.ttl1.sample_input()
            #self.mutate_dataset(self.dataset_name, 0, self.ttl1.sample_get())
            #self.core.break_realtime()
            #self.ttl4.off()
            #self.core.break_realtime()

    def analyze(self):
        th1 = self.get_dataset(self.dataset_name)
        print(th1)
        th2 = np.where(th1 > 0)
        print(th2)