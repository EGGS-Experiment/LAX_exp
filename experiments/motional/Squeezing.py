import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout,
                                         SqueezeConfigurable)


class Squeezing(LAXExperiment, Experiment):
    """
    Experiment: Squeezing

    Create a squeezed motional state by modulating the trap RF via a DDS signal
    and readout via either sideband comparison or sideband rabi flopping.
    """
    name = 'Squeezing'
    kernel_invariants = {
        'freq_sideband_readout_ftw_list', 'time_readout_mu_list',
        'freq_squeeze_ftw_list', 'phase_antisqueeze_pow_list', 'time_squeeze_mu_list', 'time_delay_mu_list',
        'config_squeeze_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=10, ndecimals=0, step=1, min=1, max=100000))

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandcool_subsequence =     SidebandCoolContinuous(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # squeezing configuration
        self.setattr_argument("freq_squeeze_khz_list",          Scannable(
                                                                    default=[
                                                                        ExplicitScan([1542.2]),
                                                                        CenterScan(1542.2, 10, 0.25, randomize=True)
                                                                    ],
                                                                    global_min=0, global_max=400000, global_step=1,
                                                                    unit="kHz", scale=1, ndecimals=5
                                                                ), group=self.name)
        self.setattr_argument("phase_antisqueeze_turns_list",   Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        RangeScan(0, 1.0, 6, randomize=True)
                                                                    ],
                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                    unit="turns", scale=1, ndecimals=3
                                                                ), group=self.name)
        self.setattr_argument("time_squeeze_us_list",           Scannable(
                                                                    default=[
                                                                        ExplicitScan([50]),
                                                                        RangeScan(8, 250, 100, randomize=True)
                                                                    ],
                                                                    global_min=2, global_max=10000000, global_step=10,
                                                                    unit="us", scale=1, ndecimals=3
                                                                ), group=self.name)
        self.setattr_argument("time_delay_us_list",             Scannable(
                                                                    default=[
                                                                        ExplicitScan([5]),
                                                                        RangeScan(5, 7, 51, randomize=True)
                                                                    ],
                                                                    global_min=2.5, global_max=10000000, global_step=5,
                                                                    unit="us", scale=1, ndecimals=3
                                                                ), group=self.name)
        self.setattr_argument("time_readout_us_list",           Scannable(
                                                                    default=[
                                                                        ExplicitScan([95]),
                                                                        RangeScan(3, 20, 100, randomize=True)
                                                                    ],
                                                                    global_min=1, global_max=100000, global_step=1,
                                                                    unit="us", scale=1, ndecimals=3
                                                                ), group=self.name)

        # set up squeezing
        self.squeeze_subsequence = SqueezeConfigurable(self)
        self.setattr_device('dds_parametric')
        self.setattr_device("qubit")
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

    def prepare_experiment(self):
        # ensure delay time is above minimum value
        if min(list(self.time_delay_us_list)) <= 1:
            raise Exception("Error: Delay time must be greater than 1 us.")

        '''READOUT PARAMETERS'''
        # get readout values
        self.freq_sideband_readout_ftw_list =   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list
        self.time_readout_mu_list =             np.array([self.core.seconds_to_mu(time_us * us)
                                                          for time_us in self.time_readout_us_list])

        '''SQUEEZING PARAMETERS'''
        self.freq_squeeze_ftw_list =        np.array([hz_to_ftw(freq_khz * kHz)
                                                      for freq_khz in self.freq_squeeze_khz_list])
        self.phase_antisqueeze_pow_list =   np.array([self.dds_parametric.turns_to_pow(phase_turns)
                                                      for phase_turns in self.phase_antisqueeze_turns_list])
        self.time_squeeze_mu_list =         np.array([self.core.seconds_to_mu(time_us * us)
                                                      for time_us in self.time_squeeze_us_list])
        # note: 2.381 is inherent system overhead; will be smaller if we stop doing stuff with urukul1_ch2
        # todo: change the inherent system overhead
        self.time_delay_mu_list =           np.array([self.core.seconds_to_mu((time_us - 2.381) * us)
                                                      for time_us in self.time_delay_us_list])

        # create an array of values for the experiment to sweep
        # (i.e. squeeze frequency, squeeze phase, squeeze time, delay time, readout FTW)
        self.config_squeeze_list = np.stack(np.meshgrid(self.freq_squeeze_ftw_list,
                                                        self.phase_antisqueeze_pow_list,
                                                        self.time_squeeze_mu_list,
                                                        self.freq_sideband_readout_ftw_list,
                                                        self.time_delay_mu_list,
                                                        self.time_readout_mu_list),
                                            -1).reshape(-1, 6)
        np.random.shuffle(self.config_squeeze_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_squeeze_list),
                7)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            for config_vals in self.config_squeeze_list:
                '''CONFIGURE'''
                # extract values from config list
                freq_squeeze_ftw =      np.int32(config_vals[0])
                phase_antisqueeze_pow = np.int32(config_vals[1])
                time_squeeze_mu =       config_vals[2]
                freq_readout_ftw =      np.int32(config_vals[3])
                time_delay_mu =         config_vals[4]
                time_readout_mu =       config_vals[5]
                self.core.break_realtime()

                # configure squeezing and qubit readout
                # note - need to return squeezing time since it is liable to change
                time_squeeze_mu = self.squeeze_subsequence.configure(freq_squeeze_ftw, phase_antisqueeze_pow, time_squeeze_mu)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                '''SQUEEZING'''
                # squeeze!
                self.squeeze_subsequence.squeeze()
                # configurable delay to simulate a pulse sequence
                delay_mu(time_delay_mu)
                # antisqueeze!
                self.squeeze_subsequence.antisqueeze()

                '''READOUT'''
                # sideband shelve
                self.sidebandreadout_subsequence.run_time(time_readout_mu)
                # read out fluorescence
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(),
                                    freq_squeeze_ftw, phase_antisqueeze_pow,
                                    time_squeeze_mu, time_delay_mu, time_readout_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    def analyze(self):
        pass
