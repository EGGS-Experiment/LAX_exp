import numpy as np
from artiq.experiment import *


class PIDTest(EnvExperiment):
    """
    Test voltage PID servo.
    """
    # todo: set relevant kernel invariants for constant values
    # todo: convert all values to appropriate voltage

    def build(self):
        # devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")

        # timing
        self.setattr_argument("time_delay_us", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))

        # PID
        self.setattr_argument("param_p", NumberValue(default=0.1, ndecimals=5, step=0.001, min=-100, max=100))
        self.setattr_argument("param_i", NumberValue(default=6.0, ndecimals=5, step=1, min=-100, max=100))
        self.setattr_argument("param_d", NumberValue(default=0, ndecimals=0, step=1, min=1, max=1000))

        # num points
        self.setattr_argument("num_points", NumberValue(default=500, ndecimals=0, step=1, min=1, max=10000))

    def prepare(self):
        """
        todo
        """
        # reassign device names
        self.adc = self.sampler0
        self.dac = self.fastino0

        # ADC
        self.channel_feedback = 0
        self.channel_setpoint = 1

        # DAC channels
        self.channel_output = 0

        # timing
        # split time_delay_mu into two
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # PID
        self.error_integral = 0
        self.err_signal = 0
        self.time_constant = self.time_delay_us * us

        # datasets
        self.set_dataset('pid_dataset', np.zeros(self.num_points), broadcast=True)
        self.setattr_dataset('pid_dataset')

        self.set_dataset('err_dataset', np.zeros(self.num_points), broadcast=True)
        self.setattr_dataset('err_dataset')

        # todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w
        # todo: may have to put r/w onto core dma

    @kernel(flags={"fast-math"})
    def run(self):
        # prepare devices
        self.core.reset()
        self.prepareDevices()
        self.core.break_realtime()

        # PID loop
        sampler_buffer = [0] * 8
        for i in range(self.num_points):
            with parallel:
                delay_mu(self.time_delay_mu)

                with sequential:
                    # get error from sampler
                    self.sampler0.sample_mu(sampler_buffer)
                    setpoint_mu = int(sampler_buffer[self.channel_setpoint] / 10) + 0x8000
                    err_val_mu = (sampler_buffer[self.channel_setpoint] - sampler_buffer[self.channel_feedback]) / 10

                    # create and record error signal
                    self.error_integral += np.int32(err_val_mu * self.time_constant)

                    with parallel:
                        # update fastino voltage
                        with sequential:
                            self.err_signal = np.int32(self.param_p * err_val_mu + self.param_i * self.error_integral) # todo: properly convert to volt_mu
                            self.core.break_realtime()
                            self.fastino0.set_dac_mu(self.channel_output, self.err_signal + setpoint_mu) # todo: properly convert to volt_mu

                        # store data
                        self.mutate_dataset('pid_dataset', i, err_val_mu)
                        self.mutate_dataset('err_dataset', i, self.err_signal)

                    # todo: consider async storing error value for display

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up sampler
        self.adc.init()
        self.core.break_realtime()
        self.adc.set_gain_mu(self.channel_feedback, 1)
        self.adc.set_gain_mu(self.channel_setpoint, 1)
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

    def analyze(self):
        """
        Print results.
        """
        print("Sampler error:")
        print(self.pid_dataset)

        print("PID Output Signal:")
        print(self.err_dataset)
