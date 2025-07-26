import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM,
    QubitRAP, Readout, RescueIon
)


class RapidAdiabaticPassage(LAXExperiment, Experiment):
    """
    Experiment: Rapid Adiabatic Passage

    Applies a chirped and pulse-shaped 729nm qubit pulse to achieve Rapid Adiabatic Passage.
    Demonstrated via rabi flopping/spectrum scanning.
    """
    name = 'Rapid Adiabatic Passage'
    kernel_invariants = {
        # hardware parameters
        'att_qubit_mu', 'ampl_pulse_readout_asf', 'att_pulse_readout_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'rap_subsequence',
        'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_readout', 'profile_729_SBC', 'profile_729_RAP',
        'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=10, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("cooling_type",   EnumerationValue(["Doppler", "SBC - Continuous"], default="Doppler"))

        # allocate relevant beam profiles
        self.profile_729_readout =  0
        self.profile_729_SBC =      1
        self.profile_729_RAP =      2

        # build SBC here so its arguments come first
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )

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
        self.setattr_argument("freq_rap_dev_khz_list",  Scannable(
                                                            default=[
                                                                ExplicitScan([100.]),
                                                                RangeScan(100, 500., 401, randomize=True),
                                                                CenterScan(250., 100., 5., randomize=True),
                                                            ],
                                                            global_min=0.1, global_max=100000, global_step=5.,
                                                            unit="kHz", scale=1, precision=6
                                                        ), group="{}.chirp".format(self.name))
        self.setattr_argument("time_rap_us_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([200.]),
                                                            RangeScan(1, 50, 200, randomize=True),
                                                            CenterScan(250., 100., 5., randomize=True),
                                                        ],
                                                        global_min=1, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="{}.chirp".format(self.name))
        self.setattr_argument("enable_cutoff",      BooleanValue(default=True), group="{}.chirp".format(self.name))
        self.setattr_argument("time_cutoff_us_list",    Scannable(
                                                            default=[
                                                                ExplicitScan([200.]),
                                                                RangeScan(1, 50, 200, randomize=True),
                                                                CenterScan(250., 100., 5., randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group="{}.chirp".format(self.name))

        # pulse parameters
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=30, precision=3, step=5, min=1, max=50, scale=1., unit="%"), group="{}.pulse".format(self.name))
        self.setattr_argument("att_qubit_db",   NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5, scale=1., unit="dB"), group="{}.pulse".format(self.name))
        # todo: actually implement the enable_chirp stuff lol
        self.setattr_argument("enable_pulseshaping",    BooleanValue(default=True), group="{}.pulse".format(self.name))
        self.setattr_argument("enable_chirp",           BooleanValue(default=True), group="{}.pulse".format(self.name))

        # read out parameters
        self.setattr_argument("enable_rabiflop_readout",    BooleanValue(default=False), group="rabiflop_readout")
        self.setattr_argument("ampl_pulse_readout_pct",     NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"), group="rabiflop_readout")
        self.setattr_argument("att_pulse_readout_db",       NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"), group="rabiflop_readout")
        self.setattr_argument("freq_pulse_readout_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.9851]),
                                                                    CenterScan(101.9851, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.9801, 101.9901, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="rabiflop_readout")
        self.setattr_argument("time_pulse_readout_us_list",    Scannable(
                                                                default=[
                                                                    ExplicitScan([122.9]),
                                                                    RangeScan(0, 1500, 100, randomize=True),
                                                                ],
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="rabiflop_readout")

        # initialize other subsequences
        pulse_shape_str = 'blackman'
        if self.enable_pulseshaping is False:
            pulse_shape_str = 'square'
        self.rap_subsequence =          QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=500,
            ampl_max_pct=self.ampl_qubit_pct, pulse_shape=pulse_shape_str
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('repump_qubit')
        self.setattr_device('pump')

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
        self.att_qubit_mu =     att_to_mu(self.att_qubit_db * dB)

        # frequency parameters
        freq_rap_center_ftw_list = np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_rap_center_mhz_list])
        freq_rap_dev_ftw_list =    np.array([hz_to_ftw(freq_khz * kHz) for freq_khz in self.freq_rap_dev_khz_list])

        # timing parameters
        time_rap_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                     for time_us in self.time_rap_us_list])
        if self.enable_cutoff:
            time_cutoff_mu_list = np.array([self.core.seconds_to_mu(time_us * us)
                                            for time_us in self.time_cutoff_us_list])
        else:
            time_cutoff_mu_list = np.array([-1])

        # rabiflop readout pulse
        self.ampl_pulse_readout_asf =   self.qubit.amplitude_to_asf(self.ampl_pulse_readout_pct / 100.)
        self.att_pulse_readout_mu =     att_to_mu(self.att_pulse_readout_db * dB)

        if self.enable_rabiflop_readout:
            freq_pulse_readout_ftw_list =   np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                      for freq_mhz in self.freq_pulse_readout_mhz_list], dtype=np.int32)
            time_pulse_readout_mu_list =    np.array([self.core.seconds_to_mu(time_us * us)
                                                      for time_us in self.time_pulse_readout_us_list], dtype=np.int64)
        else:
            freq_pulse_readout_ftw_list =   np.array([-1], dtype=np.int64)
            time_pulse_readout_mu_list =    np.array([-1], dtype=np.int64)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(freq_rap_center_ftw_list,
                                                           freq_rap_dev_ftw_list,
                                                           time_rap_mu_list,
                                                           time_cutoff_mu_list,
                                                           freq_pulse_readout_ftw_list,
                                                           time_pulse_readout_mu_list),
                                               -1).reshape(-1, 6)
        # if not sweeping cutoff times, set cutoff times same as RAP times
        if not self.enable_cutoff:
            self.config_experiment_list[:, 3] = self.config_experiment_list[:, 2]
        # ensure correct type and shuffle
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                7)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_center_ftw =   np.int32(config_vals[0])
                freq_dev_ftw =      np.int32(config_vals[1])
                time_rap_mu =       config_vals[2]
                time_cutoff_mu =    config_vals[3]
                freq_readout_ftw =  np.int32(config_vals[4])
                time_readout_mu =   config_vals[5]
                self.core.break_realtime()

                # configure RAP pulse
                self.rap_subsequence.configure(time_rap_mu, freq_center_ftw, freq_dev_ftw)
                delay_mu(50000)

                '''INITIALIZE ION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # optional: sideband cool
                if self.cooling_type == "SBC - Continuous":
                    self.sidebandcool_subsequence.run_dma()

                # run RAP pulse
                self.qubit.set_att_mu(self.att_qubit_mu)
                self.rap_subsequence.run_rap(time_cutoff_mu)

                '''READ OUT & STORE RESULTS'''
                if self.enable_rabiflop_readout:
                    self.pulse_readout(time_readout_mu, freq_readout_ftw)

                # read out fluorescence
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_center_ftw, self.readout_subsequence.fetch_count(),
                                    freq_dev_ftw, time_rap_mu, time_cutoff_mu,
                                    freq_readout_ftw, time_readout_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_readout(self, time_pulse_mu: TInt64, freq_readout_ftw: TInt32) -> TNone:
        """
        Run a readout pulse.
        Arguments:
            time_pulse_mu: length of pulse (in machine units).
            freq_readout_ftw: readout frequency (set by the double pass) in FTW.
        """
        # set up relevant beam waveforms
        with parallel:
            # ensure quench is using correct profile
            self.pump.readout()

            # set up qubit readout pulse
            with sequential:
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_pulse_readout_asf, pow_=0,
                                  profile=self.profile_729_readout, phase_mode=PHASE_MODE_CONTINUOUS)
                self.qubit.set_att_mu(self.att_pulse_readout_mu)
                self.qubit.set_profile(self.profile_729_readout)
                self.qubit.cpld.io_update.pulse_mu(8)

        # quench D-5/2 to eliminate coherence problems
        self.repump_qubit.on()
        delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
        self.repump_qubit.off()

        # run readout pulse
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()

