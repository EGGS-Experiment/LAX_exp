import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import SqueezeConfigurable, Displace
import LAX_exp.experiments.diagnostics.SidebandCooling as SidebandCooling


class SqueezeDisplace(SidebandCooling.SidebandCooling):
    """
    Experiment: Squeeze and Displace

    Sandwich a displacement between a squeeze and an antisqueeze to
    amplify any motional evolution.
    """
    name = 'SqueezeDisplace'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # squeezing configuration
        self.setattr_argument("freq_squeeze_khz",                           NumberValue(default=100, ndecimals=3, step=100, min=1, max=1000000), group='squeeze')
        self.setattr_argument("phase_antisqueeze_turns",                    NumberValue(default=0., ndecimals=3, step=0.1, min=-1., max=1.), group='squeeze')
        self.setattr_argument("time_squeeze_us",                            NumberValue(default=10., ndecimals=3, step=100, min=1, max=1000000), group='squeeze')

        # displacement configuration
        self.setattr_argument("freq_displace_khz",                          NumberValue(default=10., ndecimals=3, step=100, min=1, max=1000000), group='displace')
        self.setattr_argument("phase_displace_turns",                       NumberValue(default=10., ndecimals=3, step=100, min=1, max=1000000), group='displace')
        self.setattr_argument("time_displace_us",                           NumberValue(default=10., ndecimals=3, step=100, min=1, max=1000000), group='displace')
        # todo: integrate with rabi flopping readout and sideband readout
        self.setattr_argument("time_readout_us_list",                       Scannable(
                                                                                default=[
                                                                                    ExplicitScan([110]),
                                                                                    RangeScan(3, 20, 100, randomize=True)
                                                                                ],
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, ndecimals=3
                                                                            ), group=self.name)


        # set up squeezing
        self.squeeze_subsequence =                                          SqueezeConfigurable(self)
        self.displace_subsequence =                                         Displace(self)
        self.setattr_device('dds_modulation')

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list =                               self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list

        # convert squeezing to machine units
        self.freq_squeeze_ftw =                                             self.dds_modulation.frequency_to_ftw(self.freq_squeeze_khz * kHz)
        self.phase_antisqueeze_pow =                                        self.dds_modulation.turns_to_pow(self.phase_antisqueeze_turns)
        self.time_squeeze_mu =                                              self.core.seconds_to_mu(self.time_squeeze_us * us)

        # convert displacement to machine units
        self.freq_displace_ftw =                                            self.dds_modulation.frequency_to_ftw(self.freq_displace_khz * kHz)
        self.phase_displace_pow =                                           self.dds_modulation.turns_to_pow(self.phase_displace_turns)
        self.time_displace_mu =                                             self.core.seconds_to_mu(self.time_displace_us * us)

        # readout timing
        self.time_readout_mu_list =                                         np.array([self.core.seconds_to_mu(time_us * us)
                                                                                      for time_us in self.time_readout_us_list])

        # create experiment sweep configuration
        self.config_squeezedisplace_list =                                  np.stack(np.meshgrid(self.freq_sideband_readout_ftw_list,
                                                                                                 self.time_readout_mu_list),
                                                                                         -1).reshape(-1, 2)
        np.random.shuffle(self.config_squeezedisplace_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_squeezedisplace_list),
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # configure subsequences
        time_squeeze_mu = self.squeeze_subsequence.configure(self.freq_squeeze_ftw, self.phase_antisqueeze_pow, self.time_squeeze_mu)
        self.displace_subsequence.configure(self.freq_displace_ftw, self.phase_displace_pow, self.time_displace_mu)

        for trial_num in range(self.repetitions):
            for config_vals in self.config_squeezedisplace_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =      np.int32(config_vals[0])
                time_readout_mu =       config_vals[1]
                self.core.break_realtime()

                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()


                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()
                # squeeze!
                self.squeeze_subsequence.squeeze()


                '''MANIPULATE'''
                # displace!
                self.displace_subsequence.run()


                '''READOUT'''
                # antisqueeze!
                self.squeeze_subsequence.antisqueeze()
                # set readout waveform for qubit
                self.qubit.set_profile(0)
                self.qubit.set_att_mu(self.sidebandreadout_subsequence.att_sideband_readout_mu)
                # transfer population to D-5/2 state
                self.qubit.on()
                delay_mu(time_readout_mu)
                self.qubit.off()
                # read out fluorescence
                self.readout_subsequence.run()

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), time_readout_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

    def analyze(self):
        pass
