import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import SqueezeConfigurable
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class Squeezing(SidebandCooling.SidebandCooling):
    """
    Experiment: Squeezing

    Squeezing
    *** todo: document
    """
    name = 'Squeezing'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # squeezing configuration
        self.setattr_argument("freq_squeeze_khz_list",                      Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([2176]),
                                                                                        CenterScan(2176, 10, 0.25, randomize=True)
                                                                                    ],
                                                                                    global_min=0, global_max=100000, global_step=1,
                                                                                    unit="kHz", scale=1, ndecimals=5
                                                                                ), group=self.name)
        self.setattr_argument("phase_squeeze_turns_list",                   Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([0.]),
                                                                                        RangeScan(0, 1.0, 6, randomize=True)
                                                                                    ],
                                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                                    unit="turns", scale=1, ndecimals=3
                                                                                ), group=self.name)
        self.setattr_argument("time_squeeze_us_list",                       Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([3]),
                                                                                        RangeScan(8, 250, 38, randomize=True)
                                                                                    ],
                                                                                    global_min=0.04, global_max=10000000, global_step=10,
                                                                                    unit="us", scale=1, ndecimals=3
                                                                                ), group=self.name)
        # self.setattr_argument("enable_antisqueezing",                   BooleanValue(default=True), group=self.name)


        # set up squeezing
        self.squeeze_subsequence =                                          SqueezeConfigurable(self)
        self.setattr_device('dds_modulation')
        # tmp remove
        self.setattr_device('urukul1_ch2')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        # tmp remove

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # convert squeezing to machine units
        self.freq_squeeze_ftw_list =                                            np.array([
                                                                                    hz_to_ftw(freq_khz * kHz)
                                                                                    for freq_khz in self.freq_squeeze_khz_list
                                                                                ])
        self.phase_squeeze_pow_list =                                           np.array([
                                                                                    self.dds_modulation.turns_to_pow(phase_turns)
                                                                                    for phase_turns in self.phase_squeeze_turns_list
                                                                                ])
        self.time_squeeze_mu_list =                                             np.array([
                                                                                    self.core.seconds_to_mu(time_us * us)
                                                                                    for time_us in self.time_squeeze_us_list
                                                                                ])

        # create an array of values for the experiment to sweep
        # (i.e. squeeze frequency, squeeze phase, squeeze time, readout FTW)
        self.config_squeeze_list =                                              np.stack(np.meshgrid(self.freq_squeeze_ftw_list,
                                                                                                     self.phase_squeeze_pow_list,
                                                                                                     self.time_squeeze_mu_list,
                                                                                                     self.freq_readout_ftw_list,),
                                                                                         -1).reshape(-1, 4)
        np.random.shuffle(self.config_squeeze_list)

    @property
    def results_shape(self):
        return (self.repetitions *
                len(self.freq_squeeze_ftw_list) * len(self.phase_squeeze_pow_list) *
                len(self.time_squeeze_mu_list) * len(self.freq_readout_ftw_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()

        # record custom readout sequence
        # note: this is necessary since DMA sequences will preserve urukul attenuation register
        with self.core_dma.record('_SBC_READOUT'):
            # set readout waveform for qubit
            self.qubit.set_profile(0)
            self.qubit.set_att_mu(self.att_readout_mu)

            # transfer population to D-5/2 state
            self.rabiflop_subsequence.run()

            # read out fluorescence
            self.readout_subsequence.run()


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep experiment config: squeezing frequency and readout frequency
            for config_vals in self.config_squeeze_list:

                # extract values from config list
                freq_squeeze_ftw =  np.int32(config_vals[0])
                phase_squeeze_pow = np.int32(config_vals[1])
                time_squeeze_mu =   config_vals[2]
                freq_readout_ftw =  np.int32(config_vals[3])
                self.core.break_realtime()

                # configure squeezing and qubit readout
                time_squeeze_mu = self.squeeze_subsequence.configure(freq_squeeze_ftw, phase_squeeze_pow,
                                                                     time_squeeze_mu)
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()
                # squeeze!
                self.squeeze_subsequence.squeeze()

                # todo: configurable delay?

                # antisqueeze!
                self.squeeze_subsequence.antisqueeze()
                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(),
                                        freq_squeeze_ftw, phase_squeeze_pow,
                                        time_squeeze_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

    def analyze(self):
        pass
