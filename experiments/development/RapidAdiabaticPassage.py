import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, QubitRAP, Readout, RescueIon


class RapidAdiabaticPassage(LAXExperiment, Experiment):
    """
    Experiment: Rapid Adiabatic Passage

    Applies a chirped and pulse-shaped 729nm qubit pulse to achieve Rapid Adiabatic Passage.
    Demonstrated via rabi flopping/spectrum scanning.
    """
    name = 'Rapid Adiabatic Passage'
    kernel_invariants = {
        'freq_qubit_scan_ftw', 'ampl_qubit_asf', 'att_qubit_mu',
        'initialize_subsequence', 'rabiflop_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=10, precision=0, step=1, min=1, max=100000))

        # chirp parameters
        self.setattr_argument("freq_rap_center_mhz_list",   Scannable(
                                                                default=[
                                                                    CenterScan(101.3318, 0.01, 0.0001, randomize=True),
                                                                    ExplicitScan([101.3318]),
                                                                    RangeScan(101.2318, 101.4318, 200, randomize=True),
                                                                ],
                                                                global_min=60, global_max=200, global_step=0.01,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="{}.chirp".format(self.name))
        self.setattr_argument("freq_rap_dev_khz_list",      Scannable(
                                                                default=[
                                                                    ExplicitScan([100.]),
                                                                    RangeScan(100, 500., 401, randomize=True),
                                                                    CenterScan(250., 100., 5., randomize=True),
                                                                ],
                                                                global_min=0.1, global_max=100000, global_step=5.,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="{}.chirp".format(self.name))
        self.setattr_argument("time_rap_us_list",           Scannable(
                                                                default=[
                                                                    ExplicitScan([6.05]),
                                                                    RangeScan(1, 50, 200, randomize=True),
                                                                    CenterScan(250., 100., 5., randomize=True),
                                                                ],
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="{}.chirp".format(self.name))
        self.setattr_argument("time_cutoff_us_list",        Scannable(
                                                                default=[
                                                                    ExplicitScan([6.05]),
                                                                    RangeScan(1, 50, 200, randomize=True),
                                                                    CenterScan(250., 100., 5., randomize=True),
                                                                ],
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="{}.chirp".format(self.name))

        # pulse parameters
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=30, precision=3, step=5, min=1, max=50), group="{}.pulse".format(self.name))
        self.setattr_argument("att_qubit_db",   NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5), group="{}.pulse".format(self.name))
        self.setattr_argument("enable_pulseshaping",    BooleanValue(default=True), group="{}.pulse".format(self.name))
        self.setattr_argument("enable_chirp",           BooleanValue(default=True), group="{}.pulse".format(self.name))
        # todo: actually implement the enable_pulseshaping/chirp stuff lol

        # relevant devices
        self.setattr_device('qubit')

        # tmp remove
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        # tmp remove

        # subsequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.rabiflop_subsequence =     RabiFlop(self, time_rabiflop_us=self.time_qubit_us)
        self.rap_subsequence =          QubitRAP(self, ram_profile=0, ampl_max_pct=self.ampl_qubit_pct,
                                                        num_samples=1000)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        SANTIIZE & VERIFY INPUT
        '''
        # todo: verify freq dev values are valid
        # todo: only sweep time OR cutoff; never both

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # beam parameters
        self.ampl_qubit_asf =   self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_qubit_mu =     att_to_mu(self.att_qubit_db * dB)

        # frequency parameters
        self.freq_rap_center_ftw_list = np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_rap_center_mhz_list])
        self.freq_rap_dev_ftw_list =    np.array([hz_to_ftw(freq_khz * Hz) for freq_khz in self.freq_rap_dev_khz_list])

        # timing parameters
        self.time_rap_mu_list =     np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_rap_us_list])
        self.time_cutoff_mu_list =  np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_cutoff_us_list])

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(self.freq_rap_center_ftw_list,
                                                           self.freq_rap_dev_ftw_list,
                                                           self.time_rap_mu_list,
                                                           self.time_cutoff_mu_list),
                                               -1).reshape(-1, 4)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


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
                freq_center_ftw =   np.int32(config_vals[0])
                freq_dev_ftw =      np.int32(config_vals[1])
                time_rap_mu =       config_vals[2]
                time_cutoff_mu =    config_vals[3]
                self.core.break_realtime()

                # configure RAP pulse
                self.rap_subsequence.configure(time_rap_mu, freq_center_ftw, freq_dev_ftw)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # run RAP pulse
                self.qubit.set_att_mu(self.att_qubit_mu)
                self.rap_subsequence.run_rap(time_cutoff_mu)

                # read out
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_center_ftw, self.readout_subsequence.fetch_count(),
                                    freq_dev_ftw, time_rap_mu, time_cutoff_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

