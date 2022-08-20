import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class locking_read(EnvExperiment):
    """
    Locking Read
    Read locking data
    """

    def build(self):
        # devices
        self.setattr_device('core')

        # channels
        self.setattr_argument("channel_error", NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("channel_dac", NumberValue(default=1, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("channel_blank", NumberValue(default=2, ndecimals=0, step=1, min=0, max=7))

        # gain
        self.setattr_argument("gain_error_10dB", NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("gain_dac_10dB", NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("gain_blank_10dB", NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))

        # timing
        self.setattr_argument("time_delay_us", NumberValue(default=500, ndecimals=3, step=1, min=0.2, max=100))
        self.setattr_argument("time_total_s", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=100000))


    def prepare(self):
        # ADC
        self.adc = self.get_device("sampler0")
        self.adc_mu_to_volts_error = (10 ** (1 - self.gain_error_10dB)) / (2 ** 15)
        self.adc_mu_to_volts_dac = (10 ** (1 - self.gain_dac_10dB)) / (2 ** 15)
        self.adc_mu_to_volts_blank = (10 ** (1 - self.gain_blank_10dB)) / (2 ** 15)

        # timing
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)
        self.repetitions = np.int32(self.time_total_s / self.time_delay_us * 1e6)

        # datasets
        self.set_dataset('locking_readout', np.zeros([self.repetitions, 3]), broadcast=True)
        self.setattr_dataset('locking_readout')


    @kernel
    def run(self):
        self.core.reset()

        # set up ADC
        self.adc.init()
        self.adc.set_gain_mu(self.channel_error, self.gain_error_10dB)
        self.adc.set_gain_mu(self.channel_dac, self.gain_dac_10dB)
        self.adc.set_gain_mu(self.channel_blank, self.gain_blank_10dB)
        self.core.break_realtime()

        sampler_buffer = [0] * 8
        self.core.break_realtime()

        # sampling loop
        for i in range(self.repetitions):
            with parallel:
                delay_mu(self.time_delay_mu)
                with sequential:
                    self.adc.sample_mu(sampler_buffer)
                    self.update_dataset(i, sampler_buffer)

    @rpc(flags={"async"})
    #def update_dataset(self, i, error_mu, dac_mu, val2):
    def update_dataset(self, i, arr):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.mutate_dataset("locking_readout",
                            i,
                            [arr[self.channel_error] * self.adc_mu_to_volts_error,
                             arr[self.channel_dac] * self.adc_mu_to_volts_dac,
                             arr[self.channel_blank] * self.adc_mu_to_volts_blank])


    def analyze(self):
        print(self.locking_readout)
        for channel_num in range(len(self.locking_readout[0])):
            dataset_tmp = self.locking_readout[:, channel_num]
            print('\tch {:d}: {:.3f} +/- {:.3f} mV'.format(channel_num, np.mean(dataset_tmp) * 1000, np.std(dataset_tmp) * 1000))
        #pass
