import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class LaserPowerMeasure(EnvExperiment):
    """
    Utility: Measure Laser Power

    Read Sampler values over time.
    """
    kernel_invariants = {
        'time_delay_mu',
        'repetitions'
    }


    def build(self):
        # devices
        self.setattr_device('core')

        # sampler
        self.setattr_argument("channel_gain_dict",              PYONValue({3: 100}), group='sampler')
        self.setattr_argument("sample_rate_hz",                 NumberValue(default=100, ndecimals=3, step=1, min=1, max=5100), group='sampler')
        self.setattr_argument("time_total_s",                   NumberValue(default=1, ndecimals=0, step=1, min=1, max=100000), group='sampler')

        # photodiode
        th0 = self._HasEnvironment.__dataset_mgr
        print(th0.ddb.keys())
        # todo: frequency, resistance
        self.setattr_argument("laser_wavelength_nm",            NumberValue(default=397, ndecimals=0, step=10, min=200, max=3200), group='photodiode')
        self.setattr_argument("photodiode_termination_ohms",    NumberValue(default=10000, ndecimals=1, step=1, min=10, max=1e6), group='photodiode')
        # self.setattr_argument("photodiode_model",               EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'), group='photodiode')


        # todo: get list of photodiodes under calibrations


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

        # todo: prepare responsivity values
        # from scipy.interpolate import Akima1DInterpolator
        # ampl_calib_points =                                             self.get_dataset('calibration.temperature.asf_calibration_curve_mhz_pct')
        # ampl_calib_curve =                                              Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # datasets
        self.set_dataset('results',                     np.zeros([self.repetitions, len(self.channel_list)]))
        self.setattr_dataset('results')

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
        self.mutate_dataset("results", i, data)


    def analyze(self):
        # print out statistics of results
        print('\tResults:')
        for i in self.channel_iter:
            print('\t\tCH{:d}:\t{:.3f} +/- {:.3f} mV'.format(self.channel_list[i],
                                                             np.mean(self.results[:, i]) * 1000,
                                                             np.std(self.results[:, i]) * 1000))
