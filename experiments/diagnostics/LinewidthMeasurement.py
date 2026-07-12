import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import AbsorptionProbe, RescueIon
from sipyco import pyon
# todo: unify temperature measurement
# todo: allocate profiles


class LinewidthMeasurement(LAXExperiment, Experiment):
    """
    Experiment: Linewidth Measurement

    Measures the 397nm linewidth by doing a weak probe linescan and fitting the resulting lineshape.
    This version turns off the 866nm cooling repump to reduce heating/line-broadening.
    """
    name = 'Linewidth Measurement'
    kernel_invariants = {
        'freq_probe_scan_mhz', 'freq_probe_scan_ftw', 'ampl_probe_scan_asf',
        'time_adc_holdoff_mu', 'adc_channel_gain_mu',
        'config_linewidth_measurement_list',
        'probe_subsequence', 'rescue_subsequence'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=150, precision=0, step=1, min=1, max=100000))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz", Scannable(
            default=[
                RangeScan(81, 129, 80, randomize=True),
                ExplicitScan([120.5]),
            ],
            global_min=70, global_max=200, global_step=1,
            unit="MHz", scale=1, precision=6
        ))

        # adc (sampler) recording
        self.setattr_argument("adc_channel_num",    NumberValue(default=2, precision=0, step=1, min=0, max=7))
        self.setattr_argument("adc_channel_gain",   EnumerationValue(['1', '10', '100', '1000'], default='100'))

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')
        self.setattr_device('sampler0')

        # subsequences
        self.probe_subsequence = AbsorptionProbe(self)
        self.rescue_subsequence = RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare experimental values to reduce runtime overhead.
        """
        '''CONVERT VALUES TO MACHINE UNITS'''
        # convert probe frequency scans
        self.freq_probe_scan_mhz = np.array([freq_mhz for freq_mhz in self.freq_probe_scan_mhz])
        self.freq_probe_scan_ftw = np.array([self.pump.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz])

        # tmp remove
        # set adc holdoff time to ensure adc records when probe beam is actually on
        self.time_adc_holdoff_mu = self.core.seconds_to_mu(1000 * us)
        self.adc_channel_gain_mu = int(np.log10(int(self.adc_channel_gain)))

        '''CALIBRATE BEAM POWERS'''
        # get amplitude calibration curve from dataset manager and interpolate the points
        # interpolation is necessary to allow continuous range of frequency values
        from scipy.interpolate import Akima1DInterpolator
        # ampl_calib_points = self.get_dataset('calibration.temperature.asf_calibration_curve_mhz_pct')
        ampl_calib_points = self.get_dataset('calibration.beam_power.pump_beam.asf_calibration_curve_mhz_pct')
        ampl_calib_curve = Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # verify desired scan range is serviceable by calibration values
        min_freq_mhz = np.min(self.freq_probe_scan_mhz)
        max_freq_mhz = np.max(self.freq_probe_scan_mhz)
        if min_freq_mhz < np.min(ampl_calib_points[:, 0]):
            raise Exception("Error: lower bound of frequency scan range below valid calibration range.")
        elif max_freq_mhz > np.max(ampl_calib_points[:, 0]):
            raise Exception("Error: upper bound of frequency scan range above valid calibration range.")

        # get calibrated amplitude values
        ampl_probe_scan_pct = np.array(ampl_calib_curve(self.freq_probe_scan_mhz))
        self.ampl_probe_scan_asf = np.array([pct_to_asf(ampl_pct) for ampl_pct in ampl_probe_scan_pct])

        '''CREATE EXPERIMENT CONFIG'''
        # set up probe waveform config
        self.config_linewidth_measurement_list = np.stack(
            np.array([self.freq_probe_scan_ftw, self.ampl_probe_scan_asf])).transpose()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_probe_scan_ftw) * 2,
                4)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # set ADC channel gains
        self.sampler0.set_gain_mu(self.adc_channel_num, self.adc_channel_gain_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        buffer_sampler = [0] * 8    # create buffer to hold sampler values
        read_actual = 0
        read_control = 0

        '''MAIN LOOP'''
        for trial_num in range(self.repetitions):
            for waveform_params in self.config_linewidth_measurement_list:

                # configure probe beam waveform
                freq_probe_ftw = waveform_params[0]
                ampl_probe_asf = waveform_params[1]

                self.core.break_realtime()
                self.pump.set_mu(freq_probe_ftw, asf=ampl_probe_asf, profile=1,
                                 phase_mode=PHASE_MODE_CONTINUOUS)
                delay_mu(50000)

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

                # get stored counts from PMT & clean up loop
                self.rescue_subsequence.resuscitate()
                counts_actual = self.pmt.fetch_count()
                counts_control = self.pmt.fetch_count()

                # update datasets
                self.update_results(freq_probe_ftw, 1, counts_actual, read_actual)
                self.update_results(freq_probe_ftw, 0, counts_control, read_control)

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        """
        Process resultant spectrum and attempt to fit.

        Returns:
            res_final (np.array): background_subtracted counts
        """
        '''COLLATE DATA'''
        res_signal, res_bgr, res_final = self._process_results()

        '''FIT DATA'''
        fit_x, fit_y, textbox_str = self._fit_results(res_final)

        # format dictionary for applet plotting
        plotting_results = {'x': results_plotting_x,
                            'y': results_plotting_y,
                            'fit_x': fit_x,
                            'fit_y': fit_y,
                            'subplot_titles': f'Linewidth Measurement',
                            'subplot_x_labels': 'AOM Frequency (MHz)',
                            'subplot_y_labels': 'Signal',
                            'rid': self.scheduler.rid,
                            'textbox_strs': [textbox_str],
                            }

        self.create_matplotlib_applet(plotting_results,
                                      name=f';Linewidth Measurement',
                                      group = ['plotting', 'diagnostics'])

        return res_final

    def _process_results(self):
        """
        Process the raw results and convert to useable data

        Returns:
            A tuple containing:
                - res_signal (np.array): photon counts with 866 on
                - res_bgr (np.array): background counts with 866 off
                - res_final (np.array): final counts (i.e. signal - background)
        """
        results_tmp = np.array(self.results)[:, :3]
        # convert x-axis (frequency) from frequency tuning word (FTW) to MHz
        results_tmp[:, 0] *= 1.e3 / 0xFFFFFFFF
        # separate results into signal counts and background counts
        results_tmp = groupBy(results_tmp, column_num=1)
        _reduce_func = lambda data: np.array([np.mean(data), np.std(data)])

        # format the results into [freq_mhz, mean_counts, stdev_counts]
        res_bgr = groupBy(results_tmp[0], column_num=0, reduce_func=_reduce_func)
        res_bgr = np.concatenate((np.array([list(res_bgr.keys())]).transpose(),
                                  np.array(list(res_bgr.values()))),
                                 axis=-1)
        res_signal = groupBy(results_tmp[1], column_num=0, reduce_func=_reduce_func)
        res_signal = np.concatenate((np.array([list(res_signal.keys())]).transpose(),
                                     np.array(list(res_signal.values()))),
                                    axis=-1)
        # process results into final form
        res_final = res_signal.copy()
        res_final[:, 2] -= res_bgr[:, 2]

        # save processed results to hdf5 as a dataset
        self.set_dataset('res_signal', res_signal)
        self.set_dataset('res_bgr', res_bgr)

        return res_signal, res_bgr, res_final


    def _fit_results(self, res_final):
        """
        Fit the linewidth curve with various lineshapes

        Args:
            res_final (np.array): results

        Returns:
            A tuple containing:
                - fit_x (np.array): x values used for fitting
                - fit_y (np.array): y values used for fitting
                - textbox_str (string): texted used for textbox in applet plotting
        """
        results_plotting_x = res_final[:, 0]
        results_plotting_y = res_final[:, 1]
        fit_x = np.linspace(np.min(results_plotting_x), np.max(results_plotting_x), len(results_plotting_x) * 10)
        try:
            # fit gaussian and lorentzianprofiles
            fitter_gauss = fitGaussian()
            fitter_lorentzian = fitLorentzian()
            fit_gaussian_params, fit_gaussian_err = fitter_gauss.fit(res_final[:, :2])
            fit_lorentzian_params, fit_lorentzian_err = fitter_lorentzian.fit(res_final[:, :2])
            fit_gaussian_fwmh_mhz = np.abs(2 * (2. * fit_gaussian_params[1]) ** -0.5)
            fit_gaussian_fwmh_mhz_err = np.abs(
                fit_gaussian_fwmh_mhz * (0.5 * fit_gaussian_err[1] / fit_gaussian_params[1]))

            # save results to dataset manager for dynamic experiments
            res_dj = [fit_gaussian_params, fit_gaussian_err]
            self.set_dataset('temp.linewidthmeasurement.results', res_dj, broadcast=True, persist=False, archive=False)
            self.set_dataset('temp.linewidthmeasurement.rid', self.scheduler.rid, broadcast=True, persist=False,
                             archive=False)

            # save fitted results to hdf5 as a dataset
            self.set_dataset('fit_gaussian_params', fit_gaussian_params)
            self.set_dataset('fit_gaussian_err', fit_gaussian_err)
            self.set_dataset('fit_lorentzian_params', fit_lorentzian_params)
            self.set_dataset('fit_lorentzian_err', fit_lorentzian_err)

            '''PRINT RESULTS'''
            # print out fitted parameters
            print("\tResults - Linewidth Measurement:")
            print("\t\tGaussian Fit:")
            print("\t\t\tLinecenter:\t {:.2f} +/- {:.2f} MHz".format(fit_gaussian_params[2], fit_gaussian_err[2]))
            print("\t\t\tFWHM:\t {:.2f} +/- {:.2f} MHz".format(fit_gaussian_fwmh_mhz, fit_gaussian_fwmh_mhz_err))
            print("\t\tLorentzian Fit:")
            print("\t\t\tLinecenter:\t {:.2f} +/- {:.2f} MHz".format(fit_lorentzian_params[2], fit_lorentzian_err[2]))
            print("\t\t\tFWHM:\t {:.2f} +/- {:.2f} MHz".format(fit_lorentzian_params[1], fit_lorentzian_err[1]))
            fit_y = fitter_gauss.fit_func(fit_x, *fit_gaussian_params)

            linecenter_gaussian_mhz = fit_gaussian_params[2]
            linecenter_gaussian_mhz_err = fit_gaussian_err[2]

            linewidth_gaussian_mhz = fit_gaussian_fwmh_mhz
            linewidth_gaussian_mhz_err = fit_gaussian_fwmh_mhz_err

            textbox_str = (
                rf'Linecenter: {linecenter_gaussian_mhz:.2f)} \pm {linecenter_gaussian_mhz_err:.2f} \mathrm{MHz} \n'
                rf'Linewidth {linewidth_gaussian_mhz,:.2f)} \pm {linewidth_gaussian_mhz_err:.2f} \mathrm{MHz}')
        except Exception as e:
            print("\tUnable to Find Optimal Fit for Linewidth Measurement")
            fit_y = [None] * len(fit_x)

            linecenter_gaussian_mhz = 'N\\A'
            linecenter_gaussian_mhz_err = 'N\\A'

            linewidth_gaussian_mhz = 'N\\A'
            linewidth_gaussian_mhz_err = 'N\\A'

            textbox_str = (
                rf'Linecenter: \mathrm{N / A} \n'
                rf'Linewidth \mathrm{N / A}')

        return fit_x, fit_y, textbox_str