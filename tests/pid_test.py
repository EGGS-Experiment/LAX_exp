from artiq.experiment import *
import numpy as np


class PIDTest(EnvExperiment):
    """
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
        self.setattr_argument("param_p", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_i", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_d", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        # todo: set channels for fastino and sampler

    def prepare(self):
        """
        todo
        """
        # set list to store sampler data in
        self.sampler_buffer = [0, 0, 0, 0, 0, 0, 0]
        # todo: set integrator storage
        # todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w
        # todo: may have to put r/w onto core dma
        pass

    @kernel
    def run(self):
        self.core.reset()
        # set up sampler
        self.sampler0.init()
        self.core.break_realtime()
        self.sampler0.set_gain_mu(1, 2)
        self.sampler0.set_gain_mu(2, 2)
        self.core.break_realtime()

        # set up fastino
        self.fastino0.init()
        self.core.break_realtime()
        self.fastino0.set_continuous(0)
        self.fastino0.set_dac(0, 0)
        self.core.break_realtime()

        # PID loop
        while True:
            # read in pickoff voltage and setpoint voltage simultaneously
            self.sampler0.sample_mu(self.sampler_buffer)
            # get error value
            err_val_mu = self.sampler_buffer[0] - self.sampler_buffer[1]
            # update PID
            self.updatePID(err_val_mu)
            # todo: make sure we are calling delay correctly and not blocking
            # todo: consider async storing error value for display
            # delay until next sample time

    @kernel(flags={"fast-math"})
    def updatePID(self, err_val):
        self.error_integral += (err_val * self.time_delay_us)
        err_signal = self.param_p * err_val + self.param_i * self.error_integral
        self.fastino0.set_dac_mu(np.int32(err_signal))

    def analyze(self):
        """
        todo
        """
        # todo: try to use @rpc(flags={"async"}) or @host-only to do a non-blocking labrad
        pass