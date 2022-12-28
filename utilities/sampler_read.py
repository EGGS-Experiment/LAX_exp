import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class sampler_read(EnvExperiment):
    """
    Sampler Read
    Read Sampler values over time.
    """
    kernel_invariants = {
        'time_delay_mu',
        'repetitions'
    }

    def build(self):
        # devices
        self.setattr_device('core')

        # channels & gains
        self.setattr_argument("channel_gain_dict",              PYONValue({0: 100, 1: 100, 2:100}), tooltip='gthkim')

        # timing
        self.setattr_argument("sample_rate_hz",                 NumberValue(default=5000, ndecimals=3, step=1, min=1, max=5100))
        self.setattr_argument("time_total_s",                   NumberValue(default=1, ndecimals=0, step=1, min=1, max=100000))


    def prepare(self):
        # general
        self.channel_list =                                     list(self.channel_gain_dict.keys())
        self.channel_iter =                                     list(range(len(self.channel_gain_dict)))
        self.gain_list_mu =                                     [int(np.log10(gain_mu)) for gain_mu in self.channel_gain_dict.values()]
        self.adc_mu_to_v_list =                                 np.array([10 / (2**15 * gain_mu) for gain_mu in self.channel_gain_dict.values()])

        # ADC
        self.adc =                                              self.get_device("sampler0")

        # timing
        self.time_delay_mu =                                    self.core.seconds_to_mu(1 / self.sample_rate_hz)
        self.repetitions =                                      np.int32(self.time_total_s * self.sample_rate_hz)

        # datasets
        self.set_dataset('locking_readout',                     np.zeros([self.repetitions, len(self.channel_list)]))
        self.setattr_dataset('locking_readout')

        # save parameters
        self.set_dataset('sample_rate_hz',                      self.sample_rate_hz)
        self.set_dataset('time_total_s',                        self.time_total_s)


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # set ADC channel gains
        for i in self.channel_iter:
            self.adc.set_gain_mu(self.channel_list[i], self.gain_list_mu[i])
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
        data = np.array(volts_mu_arr)[self.channel_list] * self.adc_mu_to_v_list
        self.mutate_dataset("locking_readout", i, data)


    def analyze(self):
        # print out noise values
        print(' Results:')
        for i in self.channel_iter:
            print('\tch {:d}: {:.3f} +/- {:.3f} mV'.format(self.channel_list[i], np.mean(self.locking_readout[:, i]) * 1000, np.std(self.locking_readout[:, i]) * 1000))
        #
        # # convert from mu to volts
        # locking_readout_tmp = self.locking_readout.transpose()
        # locking_readout_tmp = np.array([
        #     adc_mu_to_volt(locking_readout_tmp[ind], self.gain_list[ind])
        #     for ind in self.channel_iter]
        # ).transpose()
        # self.locking_readout = locking_readout_tmp
