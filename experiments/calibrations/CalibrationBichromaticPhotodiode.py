from email.policy import default
from turtledemo.sorting_animate import randomize

from artiq.experiment import *

from LAX_exp.language import *
from artiq.coredevice.sampler import adc_mu_to_volt


class BichromaticCalibrationPhotodiode(LAXExperiment, Experiment):

    """
    Calibration Bichromatic Amplitude Photodiode

    Read photodiode value while varying frequency to singlepass
    """
    name = "BichromaticCalibrationPhotodiode"
    # kernel_invariants = {}

    def build_experiment(self):
        self.setattr_argument('repetitions', NumberValue(default=50, min=1, max=1000000, step=1, precision=0))
        self.setattr_argument('urukul_channel', EnumerationValue(['singlepass0', 'singlepass1', 'singlepass2']))

        # singlepass params
        self.setattr_argument('freq_singlepass_MHz_list', Scannable(default=[
                            ExplicitScan([120.339]),
                            CenterScan(120.339, 5, 0.1),
                            RangeScan(115., 125., 100)]
                            ,
                            global_min=100.,global_max=150., global_step=0.1,
                            precision=3, scale=1,  unit='MHz'), group = 'singlepass_params')

        self.setattr_argument('ampl_singlepass_pct', NumberValue(default=50., max=100., min=0.01,
                                                                           step=5, precision=2, scale=1., unit='%'), group = 'singlepass_params')

        self.setattr_argument('att_singlepass_dB', NumberValue(default=8., min=5., max=31., scale=1,
                                                              step=0.5, precision=1, unit='dB'), group = 'singlepass_params')


        self.setattr_argument("time_pulse_us",
                              NumberValue(default=100, precision=2, step=1, min=0.04, max=10000000, unit="us",
                                          scale=1.),
                              group="singlepass_params")

        # photodiode params
        self.setattr_argument('sampler_channel', NumberValue(default=0, min=0, max=7, step=1, precision=0),
                              group='sampler_params')

        self.setattr_argument('sampler_gain', EnumerationValue(['1x', '10x', '100x', '1000x']),
                              group='sampler_params')

        self.setattr_argument('sample_rate_kHz', NumberValue(default=100., min=1., max=1000, step=1, precision=1, scale=1.,
                                                             unit='kHz'),
                              group='sampler_params')

        self.setattr_device('sampler0')
        self.setattr_device('qubit')

    def prepare_experiment(self):
        self.ampl_singlepass_asf = self.qubit.amplitude_to_asf(self.ampl_singlepass_pct/100.)
        self.att_singlepass_mu = att_to_mu(self.att_singlepass_dB * dB)

        self.freq_singlepass_ftw_list = [self.qubit.frequency_to_ftw(freq_singlepass_MHz*MHz) for
                                        freq_singlepass_MHz in self.freq_singlepass_MHz_list]

        if self.urukul_channel == 'singlepass0':
            self.urukul_device = self.qubit.singlepass0
        elif self.urukul_channel == 'singlepass1':
            self.urukul_device = self.qubit.singlepass1
        elif self.urukul_channel == 'singlepass2':
            self.urukul_device = self.qubit.singlepass2
        else:
            raise ValueError('Unknown urukul channel')

        self.profile = 7

        if self.sampler_gain == '1x':
            self.sampler_gain_setting = 0
        elif self.sampler_gain == '10x':
            self.sampler_gain_setting = 1
        elif self.sampler_gain == '100x':
            self.sampler_gain_setting = 2
        elif self.sampler_gain == '1000x':
            self.sampler_gain_setting = 3
        else:
            raise ValueError('Unknown sampler_gain, must be 1x, 10x, 100x, or 1000x')

        time_sampler_delay_s = self.time_pulse_us*us / (self.time_pulse_us*us * self.sample_rate_kHz*kHz) - 3*us

        if time_sampler_delay_s <= 0:
            self.time_sampler_delay_mu = self.core.seconds_to_mu(0)
        else:
            self.time_sampler_delay_mu = self.core.seconds_to_mu(time_sampler_delay_s)



    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_singlepass_ftw_list), 2
                )

    @kernel(flags={'fast-math'})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()
        delay_mu(125000)

        # initialize sampler
        self.sampler0.init()

        self.sampler0.set_gain_mu(self.sampler_channel, self.sampler_gain_setting)

        self.urukul_device.cpld.set_profile(7)
        self.urukul_device.cpld.io_update.pulse_mu(8)

        self.urukul_device.set_att_mu(self.att_singlepass_mu)

        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

    @kernel(flags={'fast-math'})
    def run_main(self) -> TNone:
        self.core.break_realtime()
        delay_mu(125000)
        for trial_num in range(self.repetitions):
            for freq_ftw in self.freq_singlepass_ftw_list:
                self.core.break_realtime()

                data = [0.0]*8
                self.urukul_device.set_mu(ftw=freq_ftw, asf=self.ampl_singlepass_asf,
                                          profile=self.profile)
                self.sampler0.sample(data)
                delay_mu(self.time_sampler_delay_mu)

                self.update_results(freq_ftw,
                                    data[self.sampler_channel])

            self.check_termination()


    def analyze_experiment(self):

        data = self.get_dataset('results')
        freqs = data[:,0]
        voltages = adc_mu_to_volt(data[:,1], gain=self.sampler_gain_setting)
        print(voltages)
