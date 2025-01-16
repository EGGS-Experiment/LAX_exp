import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, QubitPulseShape, Readout, RescueIon
from sipyco import pyon


class LaserScan(LAXExperiment, Experiment):
    """
    Experiment: Laser Scan

    Does a 729nm laser scan; resets the ion(s) every shot.
    Supports sine-squared pulse shaping.
    """
    name = 'Laser Scan'
    kernel_invariants = {
        'freq_qubit_scan_ftw', 'ampl_qubit_asf', 'att_qubit_mu', 'time_qubit_mu',
        'initialize_subsequence', 'rabiflop_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'pulseshape_subsequence',
        'time_linetrig_holdoff_mu_list',
        'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=25, precision=0, step=1, min=1, max=100000))

        # linetrigger
        self.setattr_argument("enable_linetrigger", BooleanValue(default=False), group='linetrigger')
        self.setattr_argument("time_linetrig_holdoff_ms_list", Scannable(
            default=[
                ExplicitScan([0.1]),
                RangeScan(1, 3, 3, randomize=True),
            ],
            global_min=0.01, global_max=1000, global_step=1,
            unit="ms", scale=1, precision=3
        ), group='linetrigger')

        # scan parameters
        self.setattr_argument("freq_qubit_scan_mhz", Scannable(
            default=[
                CenterScan(101.3548, 0.01, 0.0001, randomize=True),
                ExplicitScan([101.4459]),
                RangeScan(1, 50, 200, randomize=True),
            ],
            global_min=60, global_max=200, global_step=0.01,
            unit="MHz", scale=1, precision=6
        ), group=self.name)
        self.setattr_argument("time_qubit_us", NumberValue(default=3500, precision=3, step=500, min=1, max=10000000),
                              group=self.name)
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50, precision=3, step=5, min=1, max=50),
                              group=self.name)
        self.setattr_argument("att_qubit_db", NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5),
                              group=self.name)

        self.setattr_argument("repetitions", NumberValue(default=20, precision=0, step=1, min=1, max=100000))

        # linetrigger
        self.setattr_argument("enable_linetrigger",     BooleanValue(default=True), group='linetrigger')
        self.setattr_argument("time_linetrig_holdoff_ms_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([0.1]),
                                                                    RangeScan(1, 3, 3, randomize=True),
                                                                ],
                                                                global_min=0.01, global_max=1000, global_step=1,
                                                                unit="ms", scale=1, precision=3
                                                            ), group='linetrigger')

        # scan parameters
        self.setattr_argument("freq_qubit_scan_mhz",    Scannable(
                                                            default=[
                                                                CenterScan(101.3318, 0.01, 0.0001, randomize=True),
                                                                ExplicitScan([101.3318]),
                                                                RangeScan(1, 50, 200, randomize=True),
                                                            ],
                                                            global_min=60, global_max=200, global_step=0.01,
                                                            unit="MHz", scale=1, precision=6
                                                        ), group=self.name)
        self.setattr_argument("time_qubit_us",  NumberValue(default=3500, precision=5, step=500, min=1, max=10000000), group=self.name)
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=20, precision=3, step=5, min=0.01, max=50), group=self.name)
        self.setattr_argument("att_qubit_db",   NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5), group=self.name)
        self.setattr_argument("enable_pulseshaping", BooleanValue(default=True), group=self.name)


        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('trigger_line')

        # tmp remove
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        # tmp remove

        # subsequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.rabiflop_subsequence =     RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.pulseshape_subsequence =   QubitPulseShape(self, ram_profile=0, ampl_max_pct=self.ampl_qubit_pct,
                                                        num_samples=1000)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # laser parameters
        self.freq_qubit_scan_ftw = np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz])
        self.ampl_qubit_asf = self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_qubit_mu = att_to_mu(self.att_qubit_db * dB)

        # timing
        self.time_qubit_mu = self.core.seconds_to_mu(self.time_qubit_us * us)

        '''
        CONFIGURE LINETRIGGERING
        '''
        if self.enable_linetrigger:
            self.time_linetrig_holdoff_mu_list = np.array([self.core.seconds_to_mu(time_ms * ms)
                                                           for time_ms in self.time_linetrig_holdoff_ms_list])
        else:
            self.time_linetrig_holdoff_mu_list = np.array([0])

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(self.freq_qubit_scan_ftw,
                                                           self.time_linetrig_holdoff_mu_list),
                                               -1).reshape(-1, 2)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)


    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # ensure DMA sequences use profile 0
        self.qubit.set_profile(0)
        # reduce attenuation/power of qubit beam to resolve lines
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # set up qubit pulse
        if self.enable_pulseshaping:
            self.pulseshape_subsequence.configure(self.time_qubit_mu)
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                # tmp remove
                # turn on rescue beams while waiting
                self.core.break_realtime()
                self.pump.rescue()
                self.repump_cooling.on()
                self.repump_qubit.on()
                self.pump.on()
                # tmp remove

                # extract values from config list
                freq_ftw = np.int32(config_vals[0])
                time_holdoff_mu = config_vals[1]
                self.core.break_realtime()

                # set frequency
                if self.enable_pulseshaping:
                    self.qubit.set_ftw(freq_ftw)
                else:
                    self.qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf, profile=0)
                self.core.break_realtime()

                # wait for linetrigger
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, time_holdoff_mu)

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # rabi flop & read out
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.run()
                else:
                    self.rabiflop_subsequence.run()
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_ftw, self.readout_subsequence.fetch_count(), time_holdoff_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    # ANALYSIS
    def analyze_experiment(self):
        """
        Fit data and guess potential spectral peaks.
        """
        peak_vals, results_tmp = process_laser_scan_results(self.results, self.time_qubit_us)

        # save results to hdf5 as a dataset
        self.set_dataset('spectrum_peaks', peak_vals)
        # save results to dataset manager for dynamic experiments
        self.set_dataset('temp.laserscan.results', peak_vals, broadcast=True, persist=False, archive=False)
        self.set_dataset('temp.laserscan.rid', self.scheduler.rid, broadcast=True)

        # print peaks to log for user convenience
        # ensure we don't have too many peaks before printing to log
        if len(peak_vals) < 5:
            print("\tPeaks - Laser Scan:")
            for peak_freq, peak_prob in peak_vals:
                print("\t\t{:.4f} MHz:\t{:.2f}".format(peak_freq, peak_prob))
        else:
            print("\tWarning: too many peaks detected.")

        # get results
        results_plotting = np.array(results_tmp)
        results_plotting_x, results_plotting_y = results_plotting.transpose()

        # format dictionary for applet plotting
        plotting_results = {'x': results_plotting_x,
                            'y': results_plotting_y,
                            'subplot_titles': f'Laser Scan',
                            'subplot_x_labels': 'AOM. Freq (MHz)',
                            'subplot_y_labels': 'D State Population',
                            'rid': self.scheduler.rid,
                            }


        self.set_dataset('temp.plotting.results_laserscan', pyon.encode(plotting_results), broadcast=True)

        # create applet
        self.ccb.issue("create_applet", f"Laser Scan",
                       '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_laserscan'
                       ' --num-subplots 1',
                       group="plotting.diagnostics")

        return results_tmp
