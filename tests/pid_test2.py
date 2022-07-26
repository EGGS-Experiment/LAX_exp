import numpy as np
from artiq.experiment import *

# todo: look into faster ways of writing and doing math; borrow suservo machinery for r/w


class PIDTest2(EnvExperiment):
    """
    PID Test2
    Test voltage PID servo.
    """

    kernel_invariants = {
        "time_delay_mu", "error_time_constant_us",
        "param_p", "param_i", "param_d",
        "dac_frac_to_volts_mu", "adc_mu_to_volts"
    }


    # MAIN
    def build(self):
        """
        Get experiment arguments and set devices.
        """
        # devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")
        self.setattr_device("sampler0")

        # ADC
        self.setattr_argument("adc_gain_10dB", NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("adc_channel_feedback", NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("adc_channel_setpoint", NumberValue(default=1, ndecimals=0, step=1, min=0, max=7))

        # DAC
        self.setattr_argument("dac_channel_output", NumberValue(default=0, ndecimals=0, step=1, min=0, max=31))

        # timing
        self.setattr_argument("time_delay_us", NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("record_data_interval", NumberValue(default=100, ndecimals=0, step=1, min=100, max=100000))

        # PID
        self.setattr_argument("param_p", NumberValue(default=1, ndecimals=5, step=0.1, min=-100, max=100))
        self.setattr_argument("param_i", NumberValue(default=0, ndecimals=5, step=0.1, min=-100, max=100))
        self.setattr_argument("param_d", NumberValue(default=0, ndecimals=5, step=0.1, min=-100, max=100))

    def prepare(self):
        """
        Prepare variables and datasets for the experiment.
        Convert all values to machine units to reduce overhead.
        """
        # todo: process PID parameters
        # todo: consider filter features we want
        # ADC
        self.adc = self.sampler0
        self.adc_mu_to_volts = (10 ** (1 - self.adc_gain_10dB)) / (2 ** 15)

        # DAC
        self.dac = self.fastino0
        #self.dac_frac_to_volts_mu = (2 ** 15) / 10
        self.dac_frac_to_volts_mu = (2 ** 14)

        # timing
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # PID
        self.error_signal_frac = 0.0
        self.error_integral_frac = 0.0
        self.error_time_constant_us = self.time_delay_us * us

        # datasets
        self.set_dataset('error_value_pct', [], broadcast=True)
        self.set_dataset('error_signal_pct', [], broadcast=True)
        self.set_dataset('sampler_in_val', [], broadcast=True)
        self.set_dataset('sampler_set_val', [], broadcast=True)
        self.set_dataset('dac_out_val', [], broadcast=True)
        self.setattr_dataset('error_value_pct')
        self.setattr_dataset('error_signal_pct')
        self.setattr_dataset('sampler_in_val')
        self.setattr_dataset('sampler_set_val')
        self.setattr_dataset('dac_out_val')

    @kernel(flags={"fast-math"})
    def run(self):
        # prepare devices
        self.core.reset()
        self.prepareDevices()
        self.core.break_realtime()

        # prepare for PID loop
        sampler_buffer = [0] * 8
        #counter = 0
        err_val_frac = 0.0
        #self.prev_val = 0

        # PID loop
        while True:
            with parallel:
                delay_mu(self.time_delay_mu)

                with sequential:
                    # get signals from sampler and calculate error
                    self.sampler0.sample_mu(sampler_buffer)
                    err_val_frac = (sampler_buffer[self.adc_channel_feedback] / sampler_buffer[self.adc_channel_setpoint] - 1)

                    # calculate error and PID
                    self.error_integral_frac += err_val_frac * self.error_time_constant_us
                    self.error_signal_frac = self.param_p * err_val_frac + self.param_i * self.error_integral_frac

                    # update fastino voltage
                    self.core.break_realtime()
                    self.fastino0.set_dac_mu(self.dac_channel_output, np.int32(self.error_signal_frac * self.dac_frac_to_volts_mu) + 0x8000)

            # store data
            # if (counter % self.record_data_interval) == 0:
            #     self.recordValues(
            #         err_val_frac, self.error_signal_frac,
            #         sampler_buffer[self.adc_channel_feedback], sampler_buffer[self.adc_channel_setpoint],
            #         np.int32(self.error_signal_frac * self.dac_frac_to_volts_mu)
            #     )
            #     counter = 0
            # counter += 1

    def analyze(self):
        """
        Print results.
        """
        # todo: averaged statistics
        pass


    # HELPERS
    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up sampler
        self.adc.init()
        self.core.break_realtime()
        self.adc.set_gain_mu(self.adc_channel_feedback, self.adc_gain_10dB)
        self.adc.set_gain_mu(self.adc_channel_setpoint, self.adc_gain_10dB)
        self.core.break_realtime()

        # set up fastino
        self.dac.init()
        self.core.break_realtime()
        #self.fastino0.set_continuous(0)
        #self.dac.write(0x25, self.dac_channel_output)
        self.dac.set_dac(self.dac_channel_output, 0)
        self.core.break_realtime()
        # todo: make setpoint an argument
        self.dac.set_dac(1, 1)
        self.core.break_realtime()

    @rpc(flags={"async"})
    def recordValues(self, error, output, sampler0_val, sampler1_val, dac_val_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset("error_value_pct", error * 100)
        self.append_to_dataset("error_signal_pct", output * 100)
        self.append_to_dataset("sampler_in_val", sampler0_val)
        self.append_to_dataset("sampler_set_val", sampler1_val)
        self.append_to_dataset("dac_out_val", dac_val_mu * (10000 / (2 ** 15)))
