import labrad
import numpy as np
from os import environ
from time import sleep
from datetime import datetime
from artiq.experiment import *


class LinewidthMeasurementCalibration(EnvExperiment):
    """
    Calibration: Linewidth Measurement

    Get amplitude scaling factors to compensate for frequency dependence.
    """
    # kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # search parameters
        self.setattr_argument("target_voltage_mv",                          NumberValue(default=40, ndecimals=3, step=1, min=-10000, max=10000))
        self.setattr_argument("target_tolerance_mv",                        NumberValue(default=1, ndecimals=3, step=1, min=0, max=1000))

        # DDS setup
        self.setattr_argument("dds_freq_mhz_list",                          Scannable(
                                                                                default=RangeScan(95, 125, 121, randomize=True),
                                                                                global_min=70, global_max=140, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=5
                                                                            ), group='DDS')
        self.setattr_argument("dds_channel_num",                            NumberValue(default=1, ndecimals=0, step=1, min=0, max=3), group='DDS')
        self.setattr_argument("dds_ampl_min_pct",                           NumberValue(default=5, ndecimals=2, step=1, min=1, max=40), group='DDS')
        self.setattr_argument("dds_ampl_max_pct",                           NumberValue(default=50, ndecimals=2, step=1, min=5, max=50), group='DDS')
        self.setattr_argument("dds_attenuation_db",                         NumberValue(default=14, ndecimals=1, step=0.5, min=14, max=31.5), group='DDS')

        # sampler setup
        self.setattr_argument("adc_channel_num",                            NumberValue(default=2, ndecimals=0, step=1, min=0, max=7), group='ADC')
        self.setattr_argument("adc_gain_num",                               EnumerationValue(['1', '10', '100', '1000'], default='100'), group='ADC')
        self.setattr_argument("adc_sample_num",                             NumberValue(default=100, ndecimals=0, step=1, min=1, max=5000), group='ADC')


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # todo: checks for ampl range, target val

        # tmp remove
        self.time_start =                                                   datetime.timestamp(datetime.now())
        # tmp remove

        # get devices
        self.adc =                                                          self.get_device("sampler0")
        self.dds =                                                          self.get_device("urukul1_ch{:d}".format(self.dds_channel_num))

        # convert DDS values
        self.dds_freq_mhz_list =                                            list(self.dds_freq_mhz_list)
        self.dds_freq_ftw_list =                                            np.array([self.dds.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.dds_freq_mhz_list])
        self.dds_ampl_min_asf =                                             self.dds.amplitude_to_asf(self.dds_ampl_min_pct / 100)
        self.dds_ampl_max_asf =                                             self.dds.amplitude_to_asf(self.dds_ampl_max_pct / 100)

        # holdoff for photodiode to respond to DDS (should be roughly 10x the photodiode rise time)
        self.dds_response_holdoff_mu =                                      self.core.seconds_to_mu(250 * us)

        # convert ADC values
        self.adc_gain_num =                                                 int(self.adc_gain_num)
        self.adc_v_to_mu =                                                  (2**15 * self.adc_gain_num) / 10
        self.adc_mu_to_v =                                                  10 / (2**15 * self.adc_gain_num)
        self.adc_gain_mu =                                                  np.int32(np.log10(self.adc_gain_num))
        self.adc_poll_delay_mu =                                            self.core.seconds_to_mu(200 * us)

        # convert target voltage values into ADC machine units
        self.target_voltage_mu =                                            np.int32((self.target_voltage_mv / 1000) * self.adc_v_to_mu)
        self.target_tolerance_mu =                                          np.int32((self.target_tolerance_mv / 1000) * self.adc_v_to_mu)

        # todo: add values to kernel invariants

        # set up local datasets
        self._iter_dataset = 0
        self.set_dataset("results",                                         np.zeros((len(self.dds_freq_mhz_list), 2)))
        self.setattr_dataset("results")


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()


        # MAIN SEQUENCE
        # scan DDS frequencies
        for freq_ftw in self.dds_freq_ftw_list:

            # tmp remove
            self.core.break_realtime()
            # print('\t\t\tNEW RUN: ', freq_ftw * 1e3 / 0xFFFFFFFF)
            self.core.break_realtime()
            # delay_mu(250000)
            # tmp remove

            # do a recursion search
            ampl_calib_asf = self._recursion_search(freq_ftw, self.dds_ampl_min_asf, self.dds_ampl_max_asf)
            self.core.break_realtime()

            # set frequency
            with parallel:
                self.core.break_realtime()
                self.update_dataset(freq_ftw, ampl_calib_asf)


    @kernel(flags={"fast-math"})
    def _recursion_search(self, freq_ftw, ampl_min_asf, ampl_max_asf):
        """
        Use recursion to conduct a binary search.
        """
        # get mid-range value
        ampl_tmp_asf = np.int32((ampl_min_asf + ampl_max_asf) / 2)

        # set dds amplitude
        self.dds.set_mu(freq_ftw, asf=ampl_tmp_asf)
        delay_mu(self.dds_response_holdoff_mu)

        # get voltage value
        volt_mu = self._adc_read()

        # return target value
        if abs(volt_mu - self.target_voltage_mu) <= self.target_tolerance_mu:
            return ampl_tmp_asf
        elif volt_mu < self.target_voltage_mu:
            return self._recursion_search(freq_ftw, ampl_tmp_asf, ampl_max_asf)
        else:
            return self._recursion_search(freq_ftw, ampl_min_asf, ampl_tmp_asf)

        return ampl_tmp_asf

    @kernel(flags={"fast-math"})
    def _adc_read(self):
        """
        Read ADC inputs and return an averaged result.
        """
        # create holding values
        sampler_running_avg_mu =    np.int32(0)
        sampler_buffer_mu =         [0] * 8
        adc_input_tmp_mu =          np.int32(0)

        # sampling loop
        for i in range(self.adc_sample_num):

            # read samples into buffer and get channel reading
            self.adc.sample_mu(sampler_buffer_mu)
            adc_input_tmp_mu = sampler_buffer_mu[self.adc_channel_num]

            # calculate simple moving average instead of storing values for efficiency
            sampler_running_avg_mu = np.int32((sampler_running_avg_mu * i + adc_input_tmp_mu) / (i + 1))
            delay_mu(self.adc_poll_delay_mu)

        # return average
        return sampler_running_avg_mu

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the calibration.
        """
        # set up DDS
        self.dds.set_att(self.dds_attenuation_db * dB)
        self.dds.cfg_sw(True)

        # set up ADC
        self.adc.set_gain_mu(self.adc_channel_num, self.adc_gain_mu)
        self.core.break_realtime()

    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, ampl_asf):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # save data to dataset
        self.mutate_dataset('results', self._iter_dataset, np.array([freq_ftw, ampl_asf]))

        # update dataset iterator
        self._iter_dataset += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # set calibration timestamp
        calib_timestamp = datetime.timestamp(datetime.now())
        self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

        # sort results
        _indices_sorted = np.argsort(self.results, axis=0)[:, 0]
        results_tmp = self.results[_indices_sorted].transpose()

        # convert results to desired output form
        calib_freq_ftw, calib_ampl_asf = results_tmp
        calib_freq_mhz = np.array(self.dds.ftw_to_frequency(calib_freq_ftw)) / MHz
        calib_ampl_frac = np.array(self.dds.asf_to_amplitude(calib_ampl_asf)) * 100
        calib_final = np.array([calib_freq_mhz, calib_ampl_frac]).transpose()

        # print run time
        time_stop = datetime.timestamp(datetime.now())
        print()
        print('\t\t\tDONE')
        print('\t\t\t\tTOTAL RUN TIME: {:.2f}'.format(time_stop-self.time_start))
        print()

        # add calibration values to dataset manager
        self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', calib_final, broadcast=True, persist=True)


        ### LABRAD UPLOAD ###
        # upload data to labrad for visualization in RealSimpleGrapher
        try:
            # create connections to labrad servers
            cxn =                   labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
            dv =                    cxn.data_vault
            cr =                    cxn.context()

            # create labrad dataset title
            date =                  datetime.now()
            dataset_title_tmp =     'DDS Amplitude Calibration'
            trunk =                 '{0:s}_{1:02d}:{2:02d}'.format(dataset_title_tmp, date.hour, date.minute)
            trunk_tmp =             ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk]

            # create labrad dataset
            dv.cd(trunk_tmp, True, context=cr)
            dv.new(
                dataset_title_tmp,
                [('DDS Frequency', 'MHz')],
                [('DDS Amplitude', 'Amplitude', 'pct')],
                context=cr
            )

            # upload data to labrad's Data Vault
            dv.add(calib_final, context=cr)
            print("\tLabRAD upload successful.")

        except Exception as e:
            print(e)
            print("Warning: unable to upload data to labrad.")
