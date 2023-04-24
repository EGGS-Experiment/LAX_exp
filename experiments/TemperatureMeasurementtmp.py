import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import AbsorptionProbe
from LAX_exp.system.subsequences import AbsorptionProbe2


class TemperatureMeasurement2(LAXExperiment, Experiment):
    """
    Experiment: Temperature Measurement2

    Measure the ion temperature by doing a weak probe linescan and fitting the resulting lineshape - 2.
    """

    name = 'Temperature Measurement'

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                            NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",                    Scannable(
                                                                            default=RangeScan(91, 129, 100, randomize=True),
                                                                            global_min=85, global_max=135, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=6
                                                                        ))

        # adc (sampler) recording
        self.setattr_argument("adc_channel_num",                        NumberValue(default=2, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("adc_channel_gain",                       EnumerationValue(['1', '10', '100', '1000'], default='100'))

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')
        self.setattr_device('sampler0')

        # subsequences
        self.probe_subsequence =                                        AbsorptionProbe2(self)


    def prepare_experiment(self):
        # convert probe frequency scans
        self.freq_probe_scan_mhz =                                      np.array([freq_mhz for freq_mhz in self.freq_probe_scan_mhz])
        self.freq_probe_scan_ftw =                                      np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz])

        # get amplitude calibration curve from dataset maanger and interpolate the points
        # interpolation is necessary to allow continuous range of frequency values
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                             self.get_dataset('calibration.temperature.asf_calibration_curve_mhz_pct')
        ampl_calib_curve =                                              Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # verify scan range is servicable by calibration values
        min_freq_mhz =                                                  np.min(self.freq_probe_scan_mhz)
        max_freq_mhz =                                                  np.max(self.freq_probe_scan_mhz)
        assert (min_freq_mhz < np.max(ampl_calib_points[:, 0])) & (min_freq_mhz > np.min(ampl_calib_points[:, 0])), "Error: lower bound of frequency scan range below valid calibration range."
        assert (max_freq_mhz < np.max(ampl_calib_points[:, 0])) & (max_freq_mhz > np.min(ampl_calib_points[:, 0])), "Error: upper bound of frequency scan range above valid calibration range."

        # get calibrated amplitude values
        ampl_probe_scan_pct =                                           np.array(ampl_calib_curve(self.freq_probe_scan_mhz))
        self.ampl_probe_scan_asf =                                      np.array([pct_to_asf(ampl_pct) for ampl_pct in ampl_probe_scan_pct])

        # set up probe waveform config
        self.waveform_probe_scan =                                      np.stack(np.array([self.freq_probe_scan_ftw, self.ampl_probe_scan_asf])).transpose()

        # tmp remove
        # set adc holdoff time to ensure adc records when probe beam is actually on
        self.time_adc_holdoff_mu =                                      self.core.seconds_to_mu(1000 * us)
        self.adc_channel_gain_mu =                                      int(np.log10(int(self.adc_channel_gain)))


    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_probe_scan_ftw) * 2,
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # set ADC channel gains
        self.sampler0.set_gain_mu(self.adc_channel_num, self.adc_channel_gain_mu)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # create buffer to hold sampler values
        buffer_sampler = [0] * 8
        read_actual = 0

        # main loop
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep frequency
            for waveform_params in self.waveform_probe_scan:
                self.core.break_realtime()

                # get waveform parameters and set probe beam frequency
                freq_ftw = waveform_params[0]
                ampl_asf = waveform_params[1]
                self.pump.set_mu(freq_ftw, asf=ampl_asf, profile=1)
                self.core.break_realtime()

                # get actual data
                counts_res = self.probe_subsequence.run()

                # update datasets
                self.core.break_realtime()
                self.update_results(freq_ftw, 1, counts_res, 0)
                self.core.break_realtime()
