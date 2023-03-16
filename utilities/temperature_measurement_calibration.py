import labrad
import numpy as np
from os import environ
from time import sleep
from datetime import datetime
from artiq.experiment import *


class TemperatureMeasurementCalibration(EnvExperiment):
    """
    Temperature Measurement Calibration

    Examines the resonance for the EGGS feedthrough and calculates the appropriate amplitude scalings.
    """
    # kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # DDS setup
        self.setattr_argument("dds_channel_num",                            NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_attenuation_db",                         NumberValue(default=14, ndecimals=1, step=0.5, min=14, max=31.5))
        self.setattr_argument("freq_temperaturemeasurement_mhz",            Scannable(
                                                                                default=RangeScan(80, 95, 200, randomize=False),
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=5
                                                                            ))
        # sampler setup
        self.setattr_argument("sampler_channel_num",                        NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("sampler_gain_num",                           EnumerationValue([1, 10, 100, 1000])
        self.setattr_argument("sampler_sample_rate_hz",                     NumberValue(default=5000, ndecimals=3, step=1, min=1, max=5100))
        self.setattr_argument("sampler_sample_time_s",                      NumberValue(default=1, ndecimals=0, step=1, min=1, max=100000))

        # search parameters
        self.setattr_argument("search_target_voltage_mv",                   NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("search_target_voltage_mv",                   NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # get devices
        self.adc =                                                          self.get_device("sampler0")
        self.dds =                                                          self.get_device("urukul1_ch{:d}".format(self.dds_channel_num))

        # convert DDS frequencies
        self.freq_temperaturemeasurement_mu =                               np.array([self.dds.frequency_to_ftw(freq_mhz * MHz)
                                                                                      for freq_mhz in self.freq_temperaturemeasurement_mhz])

        # convert ADC values
        self.adc_gain_mu =                                                  int(np.log10(self.sampler_gain_num))
        self.adc_time_delay_mu =                                            self.core.seconds_to_mu(1 / self.sample_rate_hz)
        self.adc_sample_iter =                                              np.arange(int(self.sampler_sample_time_s * self.sampler_sample_rate_hz))

        # set up local datasets
        self._iter_dataset = 0
        self.set_dataset("results", np.zeros((len(self.freq_temperaturemeasurement_mu), 2)))
        self.setattr_dataset("results")


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare phaser
        self.phaser_prepare()
        self.core.break_realtime()


        # MAIN SEQUENCE
        # for trial_num in range(self.repetitions):

        # sweep eggs rf frequencies
        for freq_mu in self.freq_temperaturemeasurement_mu:

            # add extra delay
            self.core.break_realtime()
            self.dds.set_mu(freq_mu, asf=)

            # set eggs carrier via the DUC
            at_mu(self.awg_board.get_next_frame_mu())
            self.awg_eggs.set_duc_frequency((freq_mhz - 85) * MHz)
            self.awg_board.duc_stb()
            self.core.break_realtime()

            # wait given time
            self.record_power(freq_mhz)
            self.core.break_realtime()

        # disable rf output
        self.awg_eggs.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def phaser_prepare(self):
        """
        Prepare phaser for the calibration.
        """
        # todo: prepare 397
        # todo: prepare gain for sampler
        pass

    @rpc
    def record_power(self, peak_freq_mhz):
        """
        Get and record the given number of peaks, sorted by amplitude.
        """
        # todo: recursion
        pass


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # convert dbm to volts

        # scale relative to max
        calib_timestamp = datetime.timestamp(datetime.now())
        power_dataset_tmp[:, 1] /= np.max(power_dataset_tmp[:, 1])
        self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', np.array([th0, th1]).transpose(), broadcast=True, persist=True)
        self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)
        # turn dataset into numpy array for ease of use
        #self.eggs_heating = np.array(self.eggs_heating)
