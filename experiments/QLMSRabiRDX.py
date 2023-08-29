import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import TickleDDS, TickleFastDDS
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class QLMSRabiRDX(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Rabi - RDX

    Quantum Logic Mass Spectroscopy - Rabi Flopping (RDX)
    Cool the ions to the ground state of motion via sideband cooling,
    then apply a tickle to create a coherent state to be read out via RSB/BSB comparison.
    *** todo: add documentation about being fast and time sweeping
    """
    name = 'QLMSRabiRDX'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # QLMS configuration
        self.setattr_argument("freq_qlms_rabi_mhz_list",                    Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([30]),
                                                                                        CenterScan(82, 2, 0.5, randomize=True)
                                                                                    ],
                                                                                    global_min=0, global_max=400, global_step=1,
                                                                                    unit="MHz", scale=1, ndecimals=3
                                                                                ), group=self.name)
        self.setattr_argument("phase_qlms_rabi_turns_list",                 Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([0]),
                                                                                        RangeScan(0, 1.0, 6, randomize=True)
                                                                                    ],
                                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                                    unit="turns", scale=1, ndecimals=3
                                                                                ), group=self.name)
        self.setattr_argument("time_qlms_rabi_ns_list",                     Scannable(
                                                                                    default=[
                                                                                        RangeScan(8, 250, 38, randomize=True),
                                                                                        ExplicitScan([250])
                                                                                    ],
                                                                                    global_min=8, global_max=10000, global_step=8,
                                                                                    unit="ns", scale=1, ndecimals=0
                                                                                ), group=self.name)

        # subsequences
        self.tickle_subsequence =                                               TickleFastDDS(self)

        # tmp remove
        self.setattr_device('dds_modulation')

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # convert QLMS modulation to machine units
        self.freq_qlms_rabi_ftw_list =                                          np.array([
                                                                                    hz_to_ftw(freq_mhz * MHz)
                                                                                    for freq_mhz in self.freq_qlms_rabi_mhz_list
                                                                                ])
        self.phase_qlms_rabi_pow_list =                                         np.array([
                                                                                    self.dds_modulation.turns_to_pow(phase_turns)
                                                                                    for phase_turns in self.phase_qlms_rabi_turns_list
                                                                                ])
        self.time_qlms_rabi_mu_list =                                           np.array([
                                                                                    self.core.seconds_to_mu(time_ns * ns)
                                                                                    for time_ns in self.time_qlms_rabi_ns_list
                                                                                ])

        # tmp remove
        # ROUND (not truncate) values to nearest multiple of 8 and remove duplicate entries
        self.time_qlms_rabi_mu_list = np.array(list(self.time_qlms_rabi_ns_list)) + 7
        self.time_qlms_rabi_mu_list = np.int64(np.round(self.time_qlms_rabi_mu_list)) & ~0x7
        self.time_qlms_rabi_mu_list = np.unique(self.time_qlms_rabi_mu_list)
        print('\ttimes: {}'.format(self.time_qlms_rabi_mu_list))
        print('\t\tlen: {}'.format(len(self.time_qlms_rabi_mu_list)))
        # tmp remove


        # create an array of values for the experiment to sweep
        # (i.e. DDS tickle frequency, tickle time, readout FTW)
        self.config_qlms_rabi_list =                                            np.stack(np.meshgrid(self.freq_qlms_rabi_ftw_list,
                                                                                                     self.phase_qlms_rabi_pow_list,
                                                                                                     self.time_qlms_rabi_mu_list,
                                                                                                     self.freq_readout_ftw_list,),
                                                                                         -1).reshape(-1, 4)
        np.random.shuffle(self.config_qlms_rabi_list)

    @property
    def results_shape(self):
        return (self.repetitions *
                len(self.freq_qlms_rabi_ftw_list) * len(self.phase_qlms_rabi_pow_list) *
                len(self.time_qlms_rabi_mu_list) * len(self.freq_readout_ftw_list),
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

            # sweep experiment config: tickle frequency and readout frequency
            for config_vals in self.config_qlms_rabi_list:

                # extract values from config list
                freq_qlms_ftw =     np.int32(config_vals[0])
                phase_qlms_pow =    np.int32(config_vals[1])
                time_qlms_mu =      config_vals[2]
                freq_readout_ftw =  np.int32(config_vals[3])
                self.core.break_realtime()

                # configure tickle and qubit readout
                self.tickle_subsequence.configure(freq_qlms_ftw, phase_qlms_pow, time_qlms_mu)
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)

                # set readout frequency
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # QLMS fast tickle
                self.tickle_subsequence.run()

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(),
                                        freq_qlms_ftw, phase_qlms_pow, time_qlms_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

    def analyze(self):
        pass
