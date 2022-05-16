from artiq.experiment import *
import numpy as np


class PIDTest(EnvExperiment):
    """
    Test voltage PID servo.
    """

    def build(self):
        # set devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")
        # set arguments
        self.setattr_argument("time_delay_us", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_p", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_i", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_d", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        # todo: set channels for fastino and sampler

    def prepare(self):
        """
        todo
        """
        # todo: set integrator storage
        # todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w
        # todo: may have to put r/w onto core dma
        pass

    @kernel
    def run(self):
        self.core.reset()
        # todo: set up
        # todo: initialize devices
        # todo: set sampler gains
        # todo: set fastino auto update
        # loop
            # read in pickoff voltage and setpoint voltage simultaneously
            # get error value
            # add value * delay to previous integrator
            # calculate output signal
            # update fastino with error signal
            # todo: consider async storing error value
            self.fastino0.set_dac(0, 0)

    def analyze(self):
        pass
