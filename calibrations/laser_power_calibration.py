from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE

import labrad
import numpy as np
from os import environ
from datetime import datetime


class LaserPowerCalibration(EnvExperiment):
    """
    Calibration: Laser Power

    Get amplitude scaling factors to compensate for frequency dependence.
    """
    # kernel_invariants = set()

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # search parameters
        self.setattr_argument("target_voltage_mv",          NumberValue(default=40, ndecimals=3, step=1, min=-10000, max=10000))
        self.setattr_argument("target_tolerance_mv",        NumberValue(default=1, ndecimals=3, step=1, min=0, max=1000))

        # DDS setup
        self.setattr_argument("dds_freq_mhz_list",          Scannable(
                                                                    default=RangeScan(95, 125, 121, randomize=True),
                                                                    global_min=70, global_max=140, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ), group='DDS')
        self.setattr_argument("dds_channel_num",            NumberValue(default=1, ndecimals=0, step=1, min=0, max=3), group='DDS')
        self.setattr_argument("dds_ampl_min_pct",           NumberValue(default=5, ndecimals=2, step=1, min=1, max=40), group='DDS')
        self.setattr_argument("dds_ampl_max_pct",           NumberValue(default=50, ndecimals=2, step=1, min=5, max=50), group='DDS')
        self.setattr_argument("dds_attenuation_db",         NumberValue(default=14, ndecimals=1, step=0.5, min=14, max=31.5), group='DDS')

        # sampler setup
        self.setattr_argument("adc_channel_num",            NumberValue(default=2, ndecimals=0, step=1, min=0, max=7), group='ADC')
        self.setattr_argument("adc_gain_num",               EnumerationValue(['1', '10', '100', '1000'], default='100'), group='ADC')
        self.setattr_argument("adc_sample_num",             NumberValue(default=100, ndecimals=0, step=1, min=1, max=5000), group='ADC')

        # result storage
        self.setattr_argument("save_to_dataset_manager",    BooleanValue(default=True), group='dataset')
        self.setattr_argument("dataset_name",               StringValue(default='pump_beam'), group='dataset')


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        '''GENERAL'''
        # todo: checks for ampl range, target val

        # get devices
        self.adc =  self.get_device("sampler0")
        self.dds =  self.get_device("urukul2_ch{:d}".format(self.dds_channel_num))

        # add extra slack before running a new frequency
        self.time_slack_mu =            self.core.seconds_to_mu(250 * us)
        # add holdoff period after changing amplitude to allow beam & photodiode to respond
        self.dds_response_holdoff_mu =  self.core.seconds_to_mu(1000 * us)


        '''PREPARE DDS'''
        # convert DDS values
        self.dds_freq_mhz_list =    list(self.dds_freq_mhz_list)
        self.dds_freq_ftw_list =    np.array([self.dds.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.dds_freq_mhz_list])
        self.dds_ampl_min_asf =     self.dds.amplitude_to_asf(self.dds_ampl_min_pct / 100.)
        self.dds_ampl_max_asf =     self.dds.amplitude_to_asf(self.dds_ampl_max_pct / 100.)


        '''PREPARE ADC'''
        # convert ADC values
        self.adc_gain_num =         int(self.adc_gain_num)
        self.adc_v_to_mu =          (2**15 * self.adc_gain_num) / 10
        self.adc_mu_to_v =          10 / (2**15 * self.adc_gain_num)
        self.adc_gain_mu =          np.int32(np.log10(self.adc_gain_num))
        self.adc_poll_delay_mu =    self.core.seconds_to_mu(200 * us)

        # convert target voltage values into ADC machine units
        self.target_voltage_mu =    np.int32((self.target_voltage_mv / 1000) * self.adc_v_to_mu)
        self.target_tolerance_mu =  np.int32((self.target_tolerance_mv / 1000) * self.adc_v_to_mu)

        # todo: add values to kernel invariants


        '''PREPARE DATASETS'''
        # set up local datasets
        self._iter_dataset = 0
        self.set_dataset("results", np.zeros((len(self.dds_freq_mhz_list), 2)))
        self.setattr_dataset("results")

        # record start time
        self.time_start =               datetime.timestamp(datetime.now())



    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def prepareDevices(self) -> TNone:
        """
        Prepare devices for the calibration.
        """
        self.core.reset()

        # set up DDS
        self.dds.sw.off()
        self.dds.cpld.set_profile(DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        self.dds.cpld.get_all_att_mu()
        self.dds.set_att(self.dds_attenuation_db * dB)
        self.dds.sw.on()
        self.core.break_realtime()

        # set up ADC
        self.adc.set_gain_mu(self.adc_channel_num, self.adc_gain_mu)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.prepareDevices()

        # MAIN LOOP
        for freq_ftw in self.dds_freq_ftw_list:

            # add slack
            self.core.break_realtime()
            at_mu(now_mu() + self.time_slack_mu)

            # search recursively for target amplitude
            ampl_calib_asf = self._recursion_search(freq_ftw, self.dds_ampl_min_asf, self.dds_ampl_max_asf)
            self.core.break_realtime()

            # store results
            self.update_dataset(freq_ftw, ampl_calib_asf)
            self.core.break_realtime()


    """
    HELPER FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def _recursion_search(self, freq_ftw: TInt32, ampl_min_asf: TInt32, ampl_max_asf: TInt32) -> TInt32:
        """
        Use recursion to conduct a binary search.
        """
        # recursion: get mid-range value
        # ampl_tmp_asf = np.int32((ampl_min_asf + ampl_max_asf) / 2)
        # note: use bit-shift for fast divide by two & ensuring int
        ampl_center_asf = (ampl_min_asf + ampl_max_asf) >> 1

        # set dds amplitude
        self.dds.set_mu(freq_ftw, asf=ampl_center_asf, profile=DEFAULT_PROFILE)
        # wait for beam-photodiode system to respond before reading voltage
        delay_mu(self.dds_response_holdoff_mu)
        volt_mu = self._adc_read()

        # recursion: check if we meet criteria
        # todo: does abs cause RPC? try doing explicit comparison
        if abs(volt_mu - self.target_voltage_mu) <= self.target_tolerance_mu:
            return ampl_center_asf
        # recursion: check if we should go low or high
        elif volt_mu < self.target_voltage_mu:
            return self._recursion_search(freq_ftw, ampl_center_asf, ampl_max_asf)
        else:
            return self._recursion_search(freq_ftw, ampl_min_asf, ampl_center_asf)

    @kernel(flags={"fast-math"})
    def _adc_read(self) -> TInt32:
        """
        Read ADC inputs and return the averaged result.
        """
        # create holding values
        sampler_running_avg_mu =    np.int32(0)
        sampler_buffer_mu =         [0] * 8
        adc_input_tmp_mu =          np.int32(0)

        # sampling loop (running average)
        for i in range(self.adc_sample_num):
            # read samples into buffer and get channel reading
            self.adc.sample_mu(sampler_buffer_mu)
            adc_input_tmp_mu = sampler_buffer_mu[self.adc_channel_num]

            # calculate simple moving average instead of storing values for efficiency
            sampler_running_avg_mu = np.int32((sampler_running_avg_mu * i + adc_input_tmp_mu) / (i + 1))
            delay_mu(self.adc_poll_delay_mu)

        # return average
        return sampler_running_avg_mu

    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, ampl_asf) -> TNone:
        """
        Records values via rpc to minimize kernel overhead.
        """
        # save data to dataset
        self.mutate_dataset('results', self._iter_dataset, np.array([freq_ftw, ampl_asf]))
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._iter_dataset / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)

        # update dataset iterator
        self._iter_dataset += 1


    """
    ANALYZE
    """
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        '''PROCESS RESULTS'''
        # sort results
        _indices_sorted = np.argsort(self.results, axis=0)[:, 0]
        results_tmp = self.results[_indices_sorted].transpose()

        # convert results to desired output form
        calib_freq_ftw, calib_ampl_asf = results_tmp
        calib_freq_mhz =    np.array(self.dds.ftw_to_frequency(calib_freq_ftw)) / MHz
        calib_ampl_frac =   np.array(self.dds.asf_to_amplitude(calib_ampl_asf)) * 100.
        calib_final =       np.array([calib_freq_mhz, calib_ampl_frac]).transpose()

        # print run time
        time_stop = datetime.timestamp(datetime.now())
        print('\n\t\tCALIBRATION COMPLETE.')
        print('\t\t\tRUN TIME: {:.2f}\n'.format(time_stop - self.time_start))

        # add calibration values to dataset manager
        if self.save_to_dataset_manager:
            _dataset_trunk = ['calibration', 'beam_power', self.dataset_name]
            self.set_dataset('.'.join(_dataset_trunk + ['calibration_timestamp']), time_stop, broadcast=True, persist=True)
            self.set_dataset('.'.join(_dataset_trunk + ['target_voltage_mv']), self.target_voltage_mv, broadcast=True, persist=True)
            self.set_dataset('.'.join(_dataset_trunk + ['asf_calibration_curve_mhz_pct']), calib_final, broadcast=True, persist=True)


        '''UPLOAD DATA TO LABRAD'''
        # upload data to labrad for visualization in RealSimpleGrapher
        try:
            # create connections to labrad servers
            cxn =                   labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
            dv =                    cxn.data_vault
            cr =                    cxn.context()

            # create labrad dataset title
            date =                  datetime.now()
            dataset_title =         'Laser Power Calibration'
            trunk_time =            '{0:s}_{1:02d}:{2:02d}'.format(dataset_title, date.hour, date.minute)
            trunk_dataset =         ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk_time]

            # create labrad dataset
            dv.cd(trunk_dataset, True, context=cr)
            dv.new(
                dataset_title,
                [('DDS Frequency', 'MHz')],
                [('DDS Amplitude', 'Amplitude', 'pct')],
                context=cr
            )
            # add dataset parameters
            dv.add_parameter("calibration_name",    self.dataset_name, context=cr)
            dv.add_parameter("target_voltage_mv",   self.target_voltage_mv, context=cr)
            dv.add_parameter("target_tolerance_mv", self.target_tolerance_mv, context=cr)

            # upload data to labrad Data Vault
            dv.add(calib_final, context=cr)
            print("\tLabRAD upload successful.")

        except Exception as e:
            print("Warning: unable to upload data to labrad.")
            print(repr(e))
