import numpy as np
from artiq.experiment import *


class PIDTest(EnvExperiment):
    """
    PID Test
    Test voltage PID servo.
    """
    # todo: set relevant kernel invariants for constant values

    def build(self):
        # set devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")
        # set arguments
        self.setattr_argument("time_delay_us", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_p", NumberValue(default=10, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_i", NumberValue(default=1, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_d", NumberValue(default=0, ndecimals=0, step=1, min=1, max=1000))

    def prepare(self):
        """
        todo
        """
        # ADC channels
        self.channel_feedback = 0
        self.channel_setpoint = 1
        # DAC channels
        self.channel_output = 0
        # set list to store sampler data in
        self.sampler_buffer = [0] * 8
        self.error_integral = 0
        # todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w
        # todo: may have to put r/w onto core dma
        pass

    @kernel
    def run(self):
        self.core.reset()
        # prepare devices
        self.prepareDevices()
        self.core.break_realtime()
        # PID loop
        while True:
            # read in pickoff voltage and setpoint voltage simultaneously
            self.sampler0.sample_mu(self.sampler_buffer)
            # get error value
            err_val_mu = self.sampler_buffer[self.channel_feedback] - self.sampler_buffer[self.channel_setpoint]
            # update PID
            self.updatePID(err_val_mu)
            self.core.break_realtime()
            # todo: make sure we are calling delay correctly and not blocking
            # todo: consider async storing error value for display
            # todo: delay until next sample time

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up sampler
        self.sampler0.init()
        self.core.break_realtime()
        self.sampler0.set_gain_mu(self.channel_feedback, 2)
        self.sampler0.set_gain_mu(self.channel_setpoint, 2)
        self.core.break_realtime()
        # set up fastino
        self.fastino0.init()
        self.core.break_realtime()
        # self.fastino0.set_continuous(0)
        self.fastino0.write(0x25, self.channel_output)
        self.fastino0.set_dac(self.channel_output, 0)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def updatePID(self, err_val_mu):
        # update error integral
        self.error_integral += (err_val_mu * self.time_delay_us)
        # create error signal
        err_signal = self.param_p * err_val_mu + self.param_i * self.error_integral
        self.core.break_realtime()
        # update fastino voltage
        self.fastino0.set_dac_mu(self.channel_output, np.int32(err_signal))
        # todo: manage break_realtimes and delays better to be more consistent
        # todo: could use with parallel and delay

    def analyze(self):
        """
        todo
        """
        # todo: try to use @rpc(flags={"async"}) or @host-only to do a non-blocking labrad
        pass
