import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon
from analysis.processing import *


class LaserScan(LAXExperiment, Experiment):
    """
    Experiment: Laser Scan

    Gets the number of counts as a function of frequency for a fixed time.
    """
    name = 'Laser Scan'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))

        # scan parameters
        self.setattr_argument("freq_qubit_scan_mhz",                Scannable(
                                                                        default=CenterScan(103.226, 0.02, 0.0005, randomize=True),
                                                                        global_min=60, global_max=200, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=6
                                                                    ), group=self.name)
        self.setattr_argument("time_qubit_us",                      NumberValue(default=5000, ndecimals=5, step=1, min=1, max=10000000), group=self.name)
        self.setattr_argument("ampl_qubit_pct",                     NumberValue(default=50, ndecimals=3, step=10, min=1, max=50), group=self.name)
        self.setattr_argument("att_qubit_db",                       NumberValue(default=28, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)

        # relevant devices
        self.setattr_device('qubit')

        # subsequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.rabiflop_subsequence =                                 RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # convert waveform values to machine units
        self.freq_qubit_scan_ftw =                                  np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz])
        self.ampl_qubit_asf =                                       self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_qubit_mu =                                         att_to_mu(self.att_qubit_db * dB)


    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_qubit_scan_mhz),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # ensure DMA sequences use profile 0
        self.qubit.set_profile(0)
        # reduce attenuation/power of qubit beam to resolve lines
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.rabiflop_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            self.core.break_realtime()

            # sweep frequency
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # rabi flop
                self.rabiflop_subsequence.run_dma()

                # read out
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()


    # ANALYSIS
    def analyze_experiment(self):
        """
        Fit data and guess potential spectral peaks.
        """
        # # todo: move to use processFluorescence2D
        # # create data structures for processing
        # results_tmp =           np.array(self.results)
        # probability_vals =      np.zeros(len(results_tmp))
        # counts_arr =            np.array(results_tmp[:, 1])
        #
        # # convert x-axis (frequency) from frequency tuning word (FTW) to MHz
        # results_tmp[:, 0] *=    1.e3 / 0xFFFFFFFF
        #
        #
        # # calculate fluorescence detection threshold
        # threshold_list =        findThresholdScikit(results_tmp[:, 1])
        # for threshold_val in threshold_list:
        #     probability_vals[np.where(counts_arr > threshold_val)] += 1.
        # # normalize probabilities and convert from D-state probability to S-state probability
        # results_tmp[:, 1] =     1. - probability_vals / len(threshold_list)
        #
        # # process dataset into x, y, with y being averaged probability
        # results_tmp =           groupBy(results_tmp, column_num=0, reduce_func=np.mean)
        # results_tmp =           np.array([list(results_tmp.keys()), list(results_tmp.values())]).transpose()
        #
        #
        # # calculate peak criteria from data
        # # todo: somehow relate peak height to shot noise (i.e. 1/sqrt(N))
        # # todo: maybe set min peak width of at least 2 points (? not sure if good idea)
        # # _peak_height =          np.power(self.repetitions, -0.5)
        # _peak_height =          0.2
        # _peak_thresh =          0.05
        # # peak distance criteria is set as ~8 kHz between points
        # _peak_dist =            int(4e-3 / (results_tmp[1, 0] - results_tmp[0, 0]))
        #
        # # calculate peaks from data and extract values
        # from scipy.signal import find_peaks
        # peaks, props =          find_peaks(results_tmp[:, 1], height=_peak_height, distance=_peak_dist)
        # peak_vals =             results_tmp[peaks]
        #
        # # fit sinc profile to results (only in the case of one peak)
        # if len(peaks) == 1:
        #     # get index step size in frequency (mhz)
        #     step_size_mhz = np.mean(results_tmp[1:, 0] - results_tmp[:-1, 0])
        #     freq_sinc_mhz = 1. / self.time_qubit_us
        #
        #     # get points +/- 6x the 1/f time for sinc fitting
        #     num_points_sinc = round(6. * freq_sinc_mhz / step_size_mhz)
        #     index_peak_center = peaks[0]
        #     index_min = max(0, index_peak_center - num_points_sinc)
        #     index_max = min(index_peak_center + num_points_sinc, len(results_tmp))
        #     points_tmp = results_tmp[index_min: index_max]
        #
        #     # fit sinc profile and replace spectrum peak with fitted value
        #     # note: division by 2 accounts for conversion between AOM freq. and abs. freq.
        #     fit_sinc_params, _ = fitSinc(points_tmp, self.time_qubit_us / 2.)
        #     peak_vals[0, 0] = fit_sinc_params[1]

        peak_vals, results_tmp = process_laser_scan_results(self.results)
        # save results to hdf5 as a dataset
        self.set_dataset('spectrum_peaks',  peak_vals)

        # save results to dataset manager for dynamic experiments
        self.set_dataset('temp.laserscan.results', peak_vals, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.laserscan.rid', self.scheduler.rid, broadcast=True, persist=False, archive=False)

        # print peaks to log for user convenience
        # ensure we don't have too many peaks before printing to log
        if len(peak_vals) < 5:
            print("\tPeaks - Laser Scan:")
            for peak_freq, peak_prob in peak_vals:
                print("\t\t{:.4f} MHz:\t{:.2f}".format(peak_freq, peak_prob))
        else:
            print("\tWarning: too many peaks detected.")
        return results_tmp
