import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class locking_read(EnvExperiment):
    """
    Locking Read
    Read locking data
    """
    kernel_invariants = {
        'time_delay_mu',
        'repetitions'
    }

    def build(self):
        # devices
        self.setattr_device('core')

        # channels
        self.setattr_argument("channel_list",           PYONValue([0, 1, 2, 3]))
        self.setattr_argument("gain_list",              PYONValue([2, 2, 3, 0]))

        # timing
        self.setattr_argument("sample_rate_hz",         NumberValue(default=5000, ndecimals=3, step=1, min=1, max=5100))
        self.setattr_argument("time_total_s",           NumberValue(default=10, ndecimals=0, step=1, min=1, max=100000))


    def prepare(self):
        # general
        self.channel_iter =                             list(range(len(self.channel_list)))

        # ADC
        self.adc =                                      self.get_device("sampler0")

        # timing
        self.time_delay_mu =                            self.core.seconds_to_mu(1 / self.sample_rate_hz)
        self.repetitions =                              np.int32(self.time_total_s * self.sample_rate_hz)

        # datasets
        self.set_dataset('locking_readout',             np.zeros([self.repetitions, len(self.channel_list)]))
        self.setattr_dataset('locking_readout')

        # save parameters
        self.set_dataset('sample_rate_hz',              self.sample_rate_hz)
        self.set_dataset('time_total_s',                self.time_total_s)


    @kernel(flags='fast-math')
    def run(self):
        self.core.reset()

        # set up ADC
        self.adc.init()

        # set ADC channel gains
        for i in self.channel_iter:
            self.adc.set_gain_mu(self.channel_list[i], self.gain_list[i])
            self.core.break_realtime()

        # create holding buffer
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
    def update_dataset(self, i, volts_mu_arr):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.mutate_dataset("locking_readout", i, np.array(volts_mu_arr)[self.channel_list])


    def analyze(self):
        # print out noise values
        for channel_num in range(len(self.channel_list)):
            dataset_tmp = adc_mu_to_volt(self.locking_readout[:, channel_num], self.gain_list[channel_num])
            print('\tch {:d}: {:.3f} +/- {:.3f} mV'.format(channel_num, np.mean(dataset_tmp) * 1000, np.std(dataset_tmp) * 1000))

        # convert from mu to volts
        locking_readout_tmp = self.locking_readout.transpose()
        locking_readout_tmp = np.array([
            adc_mu_to_volt(locking_readout_tmp[ind], self.gain_list[ind])
            for ind in self.channel_iter]
        ).transpose()
        self.locking_readout = locking_readout_tmp
