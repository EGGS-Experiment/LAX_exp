import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import AbsorptionProbe, RescueIon

# tmp testing
from LAX_exp.analysis import *


class Spectroscopy866nm(LAXExperiment, Experiment):
    """
    Experiment: Spectroscopy 866nm

    Measures the linewidth of the D3/2 - P1/2 transition by doing a weak probe linescan.
    This version DOES NOT state prep into the D3/2 state and turns on the 397nm beam on during the probe sequence.
    """
    name = 'Spectroscopy 866nm'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                NumberValue(default=500, precision=0, step=1, min=1, max=100000))

        # probe frequency scan
        self.setattr_argument("time_probe_us",              NumberValue(default=20, precision=3, step=5, min=5, max=10000), group=self.name)
        self.setattr_argument("freq_probe_scan_mhz",        Scannable(
                                                                default=RangeScan(101, 119, 200, randomize=True),
                                                                global_min=80, global_max=140, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group=self.name)

        # configure 397nm waveform
        self.setattr_argument("freq_397nm_mhz",             NumberValue(default=124, precision=6, step=1, min=1, max=400), group='397nm')
        self.setattr_argument("ampl_397nm_pct",             NumberValue(default=15, precision=3, step=10, min=1, max=50.), group='397nm')

        # adc (sampler) recording
        self.setattr_argument("adc_channel_num",            NumberValue(default=2, precision=0, step=1, min=0, max=7), group='adc')
        self.setattr_argument("adc_channel_gain",           EnumerationValue(['1', '10', '100', '1000'], default='100'), group='adc')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')
        self.setattr_device('sampler0')

        # subsequences
        self.probe_subsequence =                            AbsorptionProbe(self)
        self.rescue_subsequence =                           RescueIon(self)

    def prepare_experiment(self):
        # convert probe values
        self.time_probe_mu =            self.core.seconds_to_mu(self.time_probe_us * us)
        self.freq_probe_scan_mhz =      np.array(list(self.freq_probe_scan_mhz))
        self.freq_probe_scan_ftw =      np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz])

        # convert 397nm values
        self.freq_397nm_ftw =       self.pump.frequency_to_ftw(self.freq_397nm_mhz * MHz)
        self.ampl_397nm_asf =       self.pump.amplitude_to_asf(self.ampl_397nm_pct / 100.)

        # get amplitude calibration curve from dataset manager and interpolate the points
        # interpolation is necessary to allow continuous range of frequency values
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =     self.get_dataset('calibration.beam_power.repump_cooling_beam.asf_calibration_curve_mhz_pct')
        ampl_calib_curve =      Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # verify scan range is serviceable by calibration values
        min_freq_mhz =  np.min(self.freq_probe_scan_mhz)
        max_freq_mhz =  np.max(self.freq_probe_scan_mhz)
        if min_freq_mhz < np.min(ampl_calib_points[:, 0]):
            raise Exception("Error: lower bound of frequency scan range below valid calibration range.")
        elif max_freq_mhz > np.max(ampl_calib_points[:, 0]):
            raise Exception("Error: upper bound of frequency scan range above valid calibration range.")

        # get calibrated amplitude values
        ampl_probe_scan_pct =           np.array(ampl_calib_curve(self.freq_probe_scan_mhz))
        self.ampl_probe_scan_asf =      np.array([pct_to_asf(ampl_pct) for ampl_pct in ampl_probe_scan_pct])

        # set up probe waveform config
        self.waveform_probe_scan =      np.stack(np.array([self.freq_probe_scan_ftw, self.ampl_probe_scan_asf])).transpose()

        # tmp remove
        # set adc holdoff time to ensure adc records when probe beam is actually on
        self.time_adc_holdoff_mu =      self.core.seconds_to_mu(1000 * us)
        self.adc_channel_gain_mu =      int(np.log10(int(self.adc_channel_gain)))


    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_probe_scan_ftw) * 2,
                4)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # set up 397nm waveform for probe sequence
        self.pump.set_mu(self.freq_397nm_ftw, asf=self.ampl_397nm_asf, profile=1)

        # set ADC channel gains
        self.sampler0.set_gain_mu(self.adc_channel_num, self.adc_channel_gain_mu)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()

        # create buffer to hold sampler values
        buffer_sampler = [0] * 8
        read_actual = 0
        read_control = 0


        '''MAIN LOOP'''
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep frequency
            for waveform_params in self.waveform_probe_scan:
                self.core.break_realtime()

                # get waveform parameters and set probe beam frequency
                freq_ftw = waveform_params[0]
                ampl_asf = waveform_params[1]
                self.repump_cooling.set_mu(freq_ftw, asf=ampl_asf, profile=1)
                self.core.break_realtime()

                # get actual data
                with parallel:
                    # probe absorption at detuning (repump on)
                    with sequential:
                        self.repump_cooling.on()
                        self.probe_subsequence.run()

                    # record beam intensity via photodiode
                    with sequential:
                        delay_mu(self.time_adc_holdoff_mu)
                        self.sampler0.sample_mu(buffer_sampler)
                        read_actual = buffer_sampler[6]
                self.core.break_realtime()

                # get control data
                with parallel:
                    # probe absorption at detuning (repump off)
                    with sequential:
                        self.repump_cooling.off()
                        self.probe_subsequence.run()

                    # record beam intensity via photodiode
                    with sequential:
                        delay_mu(self.time_adc_holdoff_mu)
                        self.sampler0.sample_mu(buffer_sampler)
                        read_control = buffer_sampler[6]

                # get counts and update datasets
                counts_actual = self.pmt.fetch_count()
                counts_control = self.pmt.fetch_count()

                # update datasets
                with parallel:
                    self.core.break_realtime()
                    self.update_results(freq_ftw, 1, counts_actual, read_actual)
                    self.update_results(freq_ftw, 0, counts_control, read_control)

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        Process resultant spectrum and attempt to fit.
        """
        # create data structures for processing
        results_tmp =           np.array(self.results)[:, :3]
        # convert x-axis (frequency) from frequency tuning word (FTW) to MHz
        results_tmp[:, 0] *=    1.e3 / 0xFFFFFFFF

        # separate results into signal counts and background counts
        results_tmp =           groupBy(results_tmp, column_num=1)
        _reduce_func =          lambda data: np.array([np.mean(data), np.std(data)])

        # format the results into [freq_mhz, mean_counts, stdev_counts]
        res_bgr =               groupBy(results_tmp[0], column_num=0, reduce_func=_reduce_func)
        res_bgr =               np.concatenate((np.array([list(res_bgr.keys())]).transpose(),
                                                np.array(list(res_bgr.values()))),
                                               axis=-1)

        res_signal =            groupBy(results_tmp[1], column_num=0, reduce_func=_reduce_func)
        res_signal =            np.concatenate((np.array([list(res_signal.keys())]).transpose(),
                                                np.array(list(res_signal.values()))),
                                               axis=-1)

        # process results into final form
        res_final =             res_signal.copy()
        res_final[:, 2] -=      res_bgr[:, 2]

        # save processed results to hdf5 as a dataset
        self.set_dataset('res_signal',  res_signal)
        self.set_dataset('res_bgr',     res_bgr)

        # # fit gaussian, lorentzian, and voigt profiles
        # # todo: fit voigt profile
        # fit_gaussian_params, fit_gaussian_err =         fitGaussian(res_final[:, :2])
        # fit_lorentzian_params, fit_lorentzian_err =     fitLorentzian(res_final[:, :2])
        # # fit_voigt_params, fit_voigt_err =               fitVoigt(res_final[:, :2])
        # fit_gaussian_fwmh_mhz =                         np.abs(2 * (2. * fit_gaussian_params[1]) ** -0.5)
        # fit_gaussian_fwmh_mhz_err =                     np.abs(fit_gaussian_fwmh_mhz * (0.5 * fit_gaussian_err[1] / fit_gaussian_params[1]))
        #
        # # save results to dataset manager for dynamic experiments
        # res_dj = [fit_gaussian_params, fit_gaussian_err]
        # self.set_dataset('temp.linewidthmeasurement.results', res_dj, broadcast=True, persist=False, archive=False)
        # self.set_dataset('temp.linewidthmeasurement.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)
        #
        # # save fitted results to hdf5 as a dataset
        # self.set_dataset('fit_gaussian_params',         fit_gaussian_params)
        # self.set_dataset('fit_gaussian_err',            fit_gaussian_err)
        #
        # self.set_dataset('fit_lorentzian_params',       fit_lorentzian_params)
        # self.set_dataset('fit_lorentzian_err',          fit_lorentzian_err)
        #
        # # self.set_dataset('fit_voigt_params',            fit_voigt_params)
        # # # self.set_dataset('fit_voigt_err',               fit_voigt_err)
        #
        # # print out fitted parameters
        # print("\tResults - Linewidth Measurement:")
        # print("\t\tGaussian Fit:")
        # print("\t\t\tLinecenter:\t {:.2f} +/- {:.2f} MHz".format(fit_gaussian_params[2], fit_gaussian_err[2]))
        # print("\t\t\tFWHM:\t {:.2f} +/- {:.2f} MHz".format(fit_gaussian_fwmh_mhz, fit_gaussian_fwmh_mhz_err))
        # print("\t\tLorentzian Fit:")
        # print("\t\t\tLinecenter:\t {:.2f} +/- {:.2f} MHz".format(fit_lorentzian_params[2], fit_lorentzian_err[2]))
        # print("\t\t\tFWHM:\t {:.2f} +/- {:.2f} MHz".format(fit_lorentzian_params[1], fit_lorentzian_err[1]))
        # # print("\t\tVoigt Fit:")
        # # print("\t\t\tLinecenter:\t {:.3f} +/- {:.3f} MHz".format(fit_gaussian_params[2], fit_gaussian_err[2]))
        # # print("\t\t\tFWHM:\t\t {:.3f} +/- {:.3f} MHz".format(fit_gaussian_fwmh_mhz, fit_gaussian_fwmh_mhz_err))
        # return res_final
