import numpy as np
from artiq.experiment import *


class PIDTest2(EnvExperiment):
    """
    Test voltage PID servo - RDX.
    """
    # todo: set relevant kernel invariants for constant values

    def build(self):
        # set devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")
        # set arguments
        self.setattr_argument("time_delay_us", NumberValue(default=100, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_p", NumberValue(default=10, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_i", NumberValue(default=1, ndecimals=0, step=1, min=1, max=1000))
        self.setattr_argument("param_d", NumberValue(default=0, ndecimals=0, step=1, min=1, max=1000))
        # num points
        self.setattr_argument("num_points", NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

    def prepare(self):
        """
        todo
        """
        # reassign device names
        self.adc = self.sampler0
        self.dac = self.fastino0
        # ADC channels
        self.channel_feedback = 0
        self.channel_setpoint = 1
        # DAC channels
        self.channel_output = 0
        # set list to store sampler data in
        #self.sampler_buffer = [0] * 8
        self.error_integral = 0
        # convert time to mu
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)
        # create dataset for storage in case
        self.dataset_name = 'pid2_dataset'
        self.set_dataset(self.dataset_name, np.zeros(self.num_points), broadcast=True)
        self.setattr_dataset(self.dataset_name)
        # todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w
        # todo: may have to put r/w onto core dma

    @kernel
    def run(self):
        self.core.reset()
        # prepare devices
        self.prepareDevices()
        self.core.break_realtime()
        # PID loop
        sampler_buffer = [0] * 8
        for i in range(self.num_points):
            # read in pickoff voltage and setpoint voltage simultaneously
            #self.sampler0.sample_mu(self.sampler_buffer)
            # store values
            #self.mutate_dataset("pid2_dataset", i, self.sampler_buffer[:2])
            # get error value
            #err_val_mu = self.sampler_buffer[self.channel_feedback] - self.sampler_buffer[self.channel_setpoint]
            # update PID
            #self.updatePID(err_val_mu)
            #self.core.break_realtime()

            with parallel:
                delay_mu(self.time_delay_mu)
                with sequential:
                    self.sampler0.sample_mu(sampler_buffer)
                    self.core.break_realtime()
                    err_val = sampler_buffer[0] - sampler_buffer[1]
                    self.mutate_dataset(self.dataset_name, i, err_val)

            # todo: make sure we are calling delay correctly and not blocking
            # todo: consider async storing error value for display
            # todo: delay until next sample time

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up sampler
        self.adc.init()
        self.core.break_realtime()
        self.adc.set_gain_mu(self.channel_feedback, 0)
        self.adc.set_gain_mu(self.channel_setpoint, 0)
        self.core.break_realtime()
        # set up fastino
        self.dac.init()
        self.core.break_realtime()
        #self.fastino0.set_continuous(0)
        #self.dac.write(0x25, self.channel_output)
        # tmp remove
        self.dac.set_dac(0, 0.1)
        self.core.break_realtime()
        self.dac.set_dac(1, 0.2)
        self.core.break_realtime()
        # tmp remove
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
        print(self.pid2_dataset)
        # todo: try to use @rpc(flags={"async"}) or @host-only to do a non-blocking labrad
        pass
