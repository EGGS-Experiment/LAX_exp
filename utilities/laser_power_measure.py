import numpy as np
from artiq.experiment import *
from artiq.coredevice.sampler import adc_mu_to_volt


class LaserPowerMeasure(EnvExperiment):
    """
    Utility: Measure Laser Power

    Read in a photodiode signal via the Sampler and convert it into a laser power.
    """
    kernel_invariants = {
        'time_delay_mu',
        'repetitions'
    }


    def build(self):
        # devices
        self.setattr_device('core')

        # sampler
        self.setattr_argument("channel_gain_dict",              PYONValue({2: 1}), group='sampler')
        self.setattr_argument("sample_rate_hz",                 NumberValue(default=100, ndecimals=3, step=1, min=1, max=5100), group='sampler')
        self.setattr_argument("time_total_s",                   NumberValue(default=1, ndecimals=0, step=1, min=1, max=100000), group='sampler')

        # photodiode
        self.setattr_argument("laser_wavelength_nm",            NumberValue(default=397, ndecimals=0, step=10, min=200, max=3200), group='photodiode')
        self.setattr_argument("photodiode_termination_ohms",    NumberValue(default=100000, ndecimals=1, step=1, min=10, max=1e7), group='photodiode')
        self.setattr_argument("photodiode_model",               StringValue(default='DET36A2'), group='photodiode')

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

        # prepare responsivity values
        try:
            photodiode_responsivity_raw = self.get_dataset('calibration.photodiode.responsivity.{:s}'.format(self.photodiode_model))
        except Exception as e:
            raise Exception("Error: Unable to find responsivity curve of photodiode model in calibrations.")

        from scipy.interpolate import Akima1DInterpolator
        photodiode_responsivity_curve =                         Akima1DInterpolator(*(photodiode_responsivity_raw.transpose()))
        self.photodiode_responsivity_value =                    photodiode_responsivity_curve(self.laser_wavelength_nm)

        # datasets
        self.set_dataset('results',                             np.zeros([self.repetitions, len(self.channel_list)]))
        self.setattr_dataset('results')

        # save parameters
        self.set_dataset('sample_rate_hz',                      self.sample_rate_hz)
        self.set_dataset('time_total_s',                        self.time_total_s)
        self.set_dataset('laser_responsivity',                  self.photodiode_responsivity_value)


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
        # todo: replace with running average measurement
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
        # extract mean voltage and stdev
        voltage_mean_v = np.mean(self.results, axis=0)
        voltage_std_v = np.std(self.results, axis=0)

        # convert voltage to power (in mW)
        power_mean_uw = voltage_mean_v / (self.photodiode_termination_ohms * self.photodiode_responsivity_value) * 1e6
        power_std_uw = voltage_std_v / (self.photodiode_termination_ohms * self.photodiode_responsivity_value) * 1e6

        # store results in dataset
        self.set_dataset('voltage_mean_v', voltage_mean_v)
        self.set_dataset('power_mean_uw', power_mean_uw)
        self.set_dataset('power_std_uw', power_std_uw)

        # print out results
        print('\tResults:')
        for i in self.channel_iter:
            # check for photodiode saturation
            if voltage_mean_v[i] >= (10. ** (-1 * self.gain_list_mu[i])):
                print('\t\tCH{:d}:\tSATURATED'.format(self.channel_list[i]))
            # otherwise, print beam power
            else:
                print('\t\tCH{:d}:\t{:.3f} +/- {:.3f} uW'.format(self.channel_list[i],
                                                                 np.mean(power_mean_uw[i]),
                                                                 np.mean(power_std_uw[i])))
