import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE
from artiq.coredevice.sampler import adc_mu_to_volt


class LaserPowerMeasure(EnvExperiment):
    """
    Utility: Laser Power Measure

    Get the power of a given beam by measuring the photodiode voltage.
    """
    kernel_invariants = {
        "dds", "dds_response_holdoff_mu", "dds_freq_ftw", "dds_ampl_asf",
        "adc", "adc_gain_mu", "adc_poll_delay_mu"
    }

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")

        # todo: get all possible DDSs

        # DDS setup
        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_name",           EnumerationValue(list(dds_device_list), default='urukul2_ch1'))
        self.setattr_argument("dds_freq_mhz",       NumberValue(default=110., precision=4, step=5, min=0.1, max=400.), group='DDS')
        self.setattr_argument("dds_ampl_pct",       NumberValue(default=15., precision=2, step=5., min=0.1, max=50.), group='DDS')
        self.setattr_argument("dds_attenuation_db", NumberValue(default=14, precision=1, step=0.5, min=14, max=31.5), group='DDS')

        # sampler setup
        self.setattr_argument("adc_channel_num",    NumberValue(default=2, precision=0, step=1, min=0, max=7), group='ADC')
        self.setattr_argument("adc_gain_num",       EnumerationValue(['1', '10', '100', '1000'], default='100'), group='ADC')
        self.setattr_argument("adc_sample_num",     NumberValue(default=1000, precision=0, step=1, min=1, max=100000), group='ADC')

        # photodiode
        self.setattr_argument("beam_wavelength_nm",             NumberValue(default=397, precision=0, step=10, min=200, max=3200), group='photodiode')
        self.setattr_argument("photodiode_termination_ohms",    NumberValue(default=100000, precision=1, step=1, min=10, max=1e7), group='photodiode')
        self.setattr_argument("photodiode_model",               StringValue(default='DET36A2'), group='photodiode')

    def _get_dds_devices(self):
        """
        Get all valid DDS (AD9910) devices from the device_db.
        """
        def is_local_dds_device(v):
            return isinstance(v, dict) and (v.get('type') == 'local') and ('class' in v) and (v.get('class') == "AD9910")

        # get only local DDS devices from device_db
        return set([k for k, v in self.get_device_db().items() if is_local_dds_device(v)])

    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        '''GENERAL'''
        # todo: checks for ampl range, target val

        # get sampler device
        self.adc = self.get_device("sampler0")

        # get arbitrary urukul channel
        try:
            self.dds = self.get_device(self.dds_name)
        except Exception as e:
            print("Error: invalid DDS channel.")
            raise e

        # add holdoff period after changing amplitude to allow beam & photodiode to respond
        self.dds_response_holdoff_mu = self.core.seconds_to_mu(5000 * us)

        '''PREPARE DDS'''
        # convert DDS values
        self.dds_freq_ftw = self.dds.frequency_to_ftw(self.dds_freq_mhz * MHz)
        self.dds_ampl_asf = self.dds.amplitude_to_asf(self.dds_ampl_pct / 100.)

        '''PREPARE ADC/PHOTODIODE'''
        # convert ADC values
        self.adc_gain_num =         int(self.adc_gain_num)
        self.adc_gain_mu =          int(np.log10(self.adc_gain_num))
        self.adc_poll_delay_mu =    self.core.seconds_to_mu(100 * us)

        # prepare photodiode responsivity values
        try:
            photodiode_responsivity_raw = self.get_dataset(
                'calibration.photodiode.responsivity.{:s}'.format(self.photodiode_model)
            )
        except Exception as e:
            raise Exception("Error: Unable to find responsivity curve of photodiode model in calibrations.")

        # interpolate responsivity curve
        from scipy.interpolate import Akima1DInterpolator
        photodiode_responsivity_curve =         Akima1DInterpolator(*(photodiode_responsivity_raw.transpose()))
        self.photodiode_responsivity_value =    photodiode_responsivity_curve(self.beam_wavelength_nm)

        '''PREPARE RESULTS'''
        self.set_dataset("results_mu", np.zeros((1, 2), np.int64))
        self.setattr_dataset("results_mu")


    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def prepareDevices(self) -> TNone:
        """
        Prepare devices before running.
        """
        self.core.reset()

        # set up DDS profile - use profile 0 as working profile
        self.dds.sw.off()
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

        # set DDS attenuator without affecting other DDSs on board
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()
        self.dds.set_att(self.dds_attenuation_db * dB)
        self.dds.sw.on()
        delay_mu(5000)

        # set up ADC
        self.adc.set_gain_mu(self.adc_channel_num, self.adc_gain_mu)
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run the experimental sequence.
        """
        self.prepareDevices()

        # set dds amplitude
        self.dds.set_mu(self.dds_freq_ftw, asf=self.dds_ampl_asf, profile=0)
        # wait for beam/photodiode loop to respond before reading voltage
        delay_mu(self.dds_response_holdoff_mu)
        self._adc_read(0)
        self.core.break_realtime()

        # clean up
        self.dds.cpld.set_profile(DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()


    """
    HELPER FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def _adc_read(self, res_idx: TInt32) -> TNone:
        """
        Read ADC inputs and return the averaged result.
        Arguments:
            res_idx (TInt32): the index to save the result as.
        Returns:
            TList(TInt32): a list of [avg_volts_mu, stdev_volts_mu].
        """
        # create holding values
        running_avg_mu =    np.int64(0)
        running_stderr_mu = np.int64(0)
        sampler_buffer_mu = [0] * 8

        # sampling loop (running average)
        for i in range(self.adc_sample_num):
            # read samples into buffer and get channel reading
            self.adc.sample_mu(sampler_buffer_mu)
            adc_result_mu = sampler_buffer_mu[self.adc_channel_num]

            # implement avg and stdev as moving average filters for efficiency
            running_avg_mu =    np.int64((running_avg_mu * i + adc_result_mu) / (i + 1))
            running_stderr_mu = np.int64((running_stderr_mu * i + adc_result_mu ** 2) / (i + 1))
            delay_mu(self.adc_poll_delay_mu)

        # process running_std from avg of squares into actual std
        running_stderr_mu = running_stderr_mu - running_avg_mu ** 2

        # store results
        self.mutate_dataset("results_mu", res_idx, np.array([running_avg_mu, running_stderr_mu]))
        self.core.break_realtime()


    """
    ANALYZE
    """
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # get results from dataset and convert from machine units to voltage
        voltage_mean_mu, voltage_stderr_mu = self.results_mu[0]
        voltage_mean_v =    adc_mu_to_volt(voltage_mean_mu, gain=self.adc_gain_mu)
        voltage_stderr_v =  adc_mu_to_volt(np.sqrt(voltage_stderr_mu / self.adc_sample_num), gain=self.adc_gain_mu)
        # convert voltage to power (in mW)
        power_mean_uw =     voltage_mean_v / (self.photodiode_termination_ohms * self.photodiode_responsivity_value) * 1e6
        power_stderr_uw =   voltage_stderr_v / (self.photodiode_termination_ohms * self.photodiode_responsivity_value) * 1e6

        # store results in dataset
        self.set_dataset('voltage_mean_v', voltage_mean_v)
        self.set_dataset('voltage_stderr_v', voltage_stderr_v)
        self.set_dataset('power_mean_uw', power_mean_uw)
        self.set_dataset('power_stderr_uw', power_stderr_uw)

        # print out results
        print('\tResults:')
        # check for photodiode saturation
        if voltage_mean_v >= (10. / self.adc_gain_num):
            print('\t\tWARNING: SATURATED')
        else:
            print(
                '\t\t{} nm:\t{:.4g} +/- {:.2g} uW'.format(
                    self.beam_wavelength_nm,
                    np.mean(power_mean_uw),
                    np.mean(power_stderr_uw)
                )
            )
