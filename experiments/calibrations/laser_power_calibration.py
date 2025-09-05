import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

import labrad
from os import environ
from datetime import datetime
# todo: implement interpolation search or some weighted search


class LaserPowerCalibration(EnvExperiment):
    """
    Calibration: Laser Power

    Get amplitude scaling factors to compensate for frequency dependence.
    """
    kernel_invariants = {
        # devices
        "adc", "dds",

        # hardware values - DDS
        "time_slack_mu", "dds_update_delay_us", "dds_update_delay_mu", "dds_freq_mhz_list",
        "dds_freq_ftw_list", "dds_ampl_min_asf", "dds_ampl_max_asf",

        # hardware values - ADC
        "adc_gain_mu", "adc_poll_delay_us", "adc_poll_delay_mu", "target_voltage_mu",
        "target_tolerance_mu",

        # other
        "time_start", '_iter_dataset', 'max_recursions',
    }

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")

        # search parameters
        self.setattr_argument("target_voltage_mv",      NumberValue(default=55.1, precision=3, step=1, min=-10000, max=10000, scale=1., unit='mV'))
        self.setattr_argument("target_tolerance_mv",    NumberValue(default=5, precision=3, step=1, min=0, max=1000, scale=1., unit='mV'))

        # DDS setup
        self.setattr_argument("dds_freq_mhz_list",  Scannable(
                                                            default=[
                                                                RangeScan(70, 140, 351, randomize=True),
                                                            ],
                                                            global_min=70, global_max=140, global_step=1,
                                                            unit="MHz", scale=1, precision=5
                                                        ), group='DDS')

        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_target",         EnumerationValue(list(dds_device_list), default='urukul0_ch1'), group='DDS')
        self.setattr_argument("dds_ampl_min_pct",   NumberValue(default=1, precision=2, step=1, min=0.1, max=50., scale=1., unit='%'), group='DDS')
        self.setattr_argument("dds_ampl_max_pct",   NumberValue(default=50, precision=2, step=1, min=0.1, max=50., scale=1., unit='%'), group='DDS')
        self.setattr_argument("dds_attenuation_db", NumberValue(default=14, precision=1, step=0.5, min=10, max=31.5, scale=1., unit='dB'), group='DDS')

        # sampler setup
        self.setattr_argument("adc_channel_num",    NumberValue(default=2, precision=0, step=1, min=0, max=7), group='ADC')
        self.setattr_argument("adc_gain_num",       EnumerationValue(['1', '10', '100', '1000'], default='100'), group='ADC')
        self.setattr_argument("adc_sample_num",     NumberValue(default=250, precision=0, step=1, min=1, max=5000), group='ADC')

        # result storage
        self.setattr_argument("save_to_dataset_manager",    BooleanValue(default=True), group='dataset')
        self.setattr_argument("dataset_name",               StringValue(default='pump_beam'), group='dataset')

        # MAGIC NUMBERS
        self.dds_update_delay_us =  1000    # slack after DDS updates for beam & photodiode to respond
        self.adc_poll_delay_us =    200     # slack between successive ADC reads
        self.max_recursions =       14

    def _get_dds_devices(self):
        """
        Get all valid DDS (AD9910) devices from the device_db.
        """
        is_local_dds_device = lambda v: (
                isinstance(v, dict) and (v.get('type') == 'local')
                and ('class' in v) and (v.get('class') == "AD9910")
        )

        # return sorted list of local DDS devices from device_db
        return sorted(set([
            k
            for k, v in self.get_device_db().items()
            if is_local_dds_device(v)
        ]))

    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        '''GENERAL'''
        # check arguments for validity
        self._prepare_argument_checks()

        # get devices
        self.adc =  self.get_device("sampler0")
        # get target DDS device
        try:
            self.dds = self.get_device(self.dds_target)
        except Exception as e:
            raise e

        # monitor error conditions
        self._error_flag = False    # data does NOT upload to dataset manager in event of error
        self._num_recursions = 0    # count recursion depth to prevent death loops


        '''PREPARE HARDWARE VALUES'''
        # timings
        self.dds_update_delay_mu =  self.core.seconds_to_mu(self.dds_update_delay_us * us) # slack after DDS updates for beam & photodiode to respond
        self.adc_poll_delay_mu =    self.core.seconds_to_mu(self.adc_poll_delay_us * us)   # slack between successive ADC reads

        # convert DDS values
        self.dds_freq_mhz_list =    list(self.dds_freq_mhz_list)
        self.dds_freq_ftw_list =    np.array([self.dds.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.dds_freq_mhz_list])
        self.dds_ampl_min_asf =     self.dds.amplitude_to_asf(self.dds_ampl_min_pct / 100.)
        self.dds_ampl_max_asf =     self.dds.amplitude_to_asf(self.dds_ampl_max_pct / 100.)

        # convert ADC values
        adc_v_to_mu =               (1<<15) * int(self.adc_gain_num) / 10
        self.adc_gain_mu =          int(np.log10(int(self.adc_gain_num)))
        self.target_voltage_mu =    np.int32((self.target_voltage_mv / 1000) * adc_v_to_mu)
        self.target_tolerance_mu =  np.int32((self.target_tolerance_mv / 1000) * adc_v_to_mu)


        '''PREPARE DATASETS'''
        self._iter_dataset = 0
        self.set_dataset("results", np.zeros((len(self.dds_freq_mhz_list), 2)))
        self.setattr_dataset("results")
        self.time_start = datetime.timestamp(datetime.now())    # record start time

    def _prepare_argument_checks(self):
        """
        Check experiment arguments for validity.
        """
        # ensure amplitude search range is valid
        if self.dds_ampl_min_pct > self.dds_ampl_max_pct:
            raise ValueError("Invalid amplitude range: [{:f}, {:f}]."
                             "dds_ampl_min_pct must be less than dds_ampl_max_pct.".format(self.dds_ampl_min_pct, self.dds_ampl_max_pct))

        # ensure target voltage is in range and won't be saturated
        if not abs(self.target_voltage_mv) <= 1e4:
            raise ValueError("Invalid target voltage: {:f}. Must be in range [-10, 10]V.".format(self.target_voltage_mv))
        elif not abs(self.target_voltage_mv) < (1e4 / int(self.adc_gain_num)):
            raise ValueError("Selected target voltage and gain will cause saturation: {:f}, {:s}.".format(self.target_voltage_mv, self.adc_gain_num))


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def prepareDevices(self) -> TNone:
        """
        Prepare devices for the calibration.
        """
        self.core.reset()

        # set up DDS profile
        self.dds.sw.off()
        self.dds.cpld.set_profile(DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)

        # set DDS attenuator without affecting other DDSs on board
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()
        self.dds.set_att(self.dds_attenuation_db * dB)
        self.dds.sw.on()

        # configure ADC
        self.adc.set_gain_mu(self.adc_channel_num, self.adc_gain_mu)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run the experimental sequence.
        """
        self.prepareDevices()

        # MAIN LOOP
        try:
            for freq_ftw in self.dds_freq_ftw_list:
                self._num_recursions = 0    # clear recursion counter

                # search recursively for target amplitude
                self.core.break_realtime()
                ampl_calib_asf = self._recursion_search(freq_ftw, self.dds_ampl_min_asf, self.dds_ampl_max_asf)

                # store results & check termination
                self.update_dataset(freq_ftw, ampl_calib_asf)
                if self.scheduler.check_termination():
                    raise TerminationRequested

        except Exception as e:
            self._error_flag = True
            print(e)
        finally:
            self.core.break_realtime()

        # todo: cleanup


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def _recursion_search(self, freq_ftw: TInt32=0x0, ampl_min_asf: TInt32=0x0, ampl_max_asf: TInt32=0x0) -> TInt32:
        """
        Use recursion to conduct a binary search for the target amplitude.
        :param freq_ftw: the frequency to set the DDS (in ftw).
        :param ampl_min_asf: the lower bound on amplitude (in asf).
        :param ampl_max_asf: the upper bound on amplitude (in asf).
        :return: the center amplitude (in asf).
        """
        # check recursion depth is OK (to prevent death spirals)
        self._num_recursions += 1
        if self._num_recursions > self.max_recursions:
            raise ValueError("Unable to find target DDS amplitude at {:f} MHz."
                             "Max recursion depth exceeded.".format(self.dds.ftw_to_frequency(freq_ftw) / MHz))

        # get mid-range value via bit-shift for fast divide by two & ensuring int
        ampl_center_asf = (ampl_min_asf + ampl_max_asf) >> 1

        # update DDS with target waveform
        self.dds.set_mu(freq_ftw, asf=ampl_center_asf,
                        profile=DEFAULT_PROFILE,
                        phase_mode=PHASE_MODE_CONTINUOUS)
        # wait for beam & photodiode to settle before reading voltage
        delay_mu(self.dds_update_delay_mu)
        volt_mu = self._adc_read()

        # recursion: check if we meet criteria
        if abs(volt_mu - self.target_voltage_mu) <= self.target_tolerance_mu:
            return ampl_center_asf
        # recursion: check if we should go low or high
        elif volt_mu < self.target_voltage_mu:
            return self._recursion_search(freq_ftw, ampl_center_asf, ampl_max_asf)
        else:
            return self._recursion_search(freq_ftw, ampl_min_asf, ampl_center_asf)

        # note: explicit return needed for artiq compiler, otherwise error
        return ampl_center_asf

    @kernel(flags={"fast-math"})
    def _adc_read(self) -> TInt32:
        """
        Read ADC inputs and return the averaged result.
        :return: the averaged ADC reading.
        """
        # create holder variables
        sampler_running_avg_mu =    np.int32(0)
        sampler_buffer_mu =         [0] * 8
        adc_input_tmp_mu =          np.int32(0)

        # sampling loop (running average)
        self.core.break_realtime()
        for i in range(self.adc_sample_num):
            # read samples into buffer and acquire channel reading
            self.adc.sample_mu(sampler_buffer_mu)
            adc_input_tmp_mu = sampler_buffer_mu[self.adc_channel_num]

            # calculate simple moving average instead of storing values for efficiency
            sampler_running_avg_mu = np.int32((sampler_running_avg_mu * i + adc_input_tmp_mu) / (i + 1))
            delay_mu(self.adc_poll_delay_mu)

        # return average
        return sampler_running_avg_mu

    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw: TInt32, ampl_asf: TInt32) -> TNone:
        """
        Records values via rpc to minimize kernel overhead.
        :param freq_ftw: the frequency of the DDS (in ftw)
        :param ampl_asf: the amplitude of the DDS (in asf).
        """
        # save data to dataset
        self.mutate_dataset('results', self._iter_dataset, np.array([freq_ftw, ampl_asf]))
        # update progress bar
        self.set_dataset('management.dynamic.completion_pct',
                         round(100. * self._iter_dataset / len(self.results), 3),
                         broadcast=True, persist=True, archive=False)
        self._iter_dataset += 1 # update dataset iterator


    '''
    ANALYZE
    '''
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
        if (not self._error_flag) and self.save_to_dataset_manager:
            _dataset_trunk = ['calibration', 'beam_power', self.dataset_name]
            self.set_dataset('.'.join(_dataset_trunk + ['calibration_timestamp']), time_stop, broadcast=True, persist=True)
            self.set_dataset('.'.join(_dataset_trunk + ['target_voltage_mv']), self.target_voltage_mv, broadcast=True, persist=True)
            self.set_dataset('.'.join(_dataset_trunk + ['asf_calibration_curve_mhz_pct']), calib_final, broadcast=True, persist=True)


        '''UPLOAD DATA TO LABRAD'''
        try:
            # create connections to labrad servers
            cxn =   labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
            dv =    cxn.data_vault
            cr =    cxn.context()

            # create labrad dataset title
            date =          datetime.now()
            dataset_title = 'Laser Power Calibration'
            trunk_time =    '{0:s}_{1:02d}:{2:02d}'.format(dataset_title, date.hour, date.minute)
            trunk_dataset = ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk_time]

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
            print("Warning - unable to upload data to labrad: {}".format(repr(e)))

