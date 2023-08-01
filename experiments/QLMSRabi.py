import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import TickleDDS
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class QLMSRabi(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Rabi

    Quantum Logic Mass Spectroscopy - Rabi Flopping
    Cool the ions to the ground state of motion via sideband cooling,
    then apply a tickle to create a coherent state to be read out via RSB/BSB comparison.
    """
    name = 'QLMSRabi'


    def build_experiment(self):
        # QLMS configuration
        self.setattr_argument("freq_qlms_rabi_khz_list",                        Scannable(
                                                                                    default=CenterScan(1101, 10, 0.25, randomize=True),
                                                                                    global_min=0, global_max=10000, global_step=1,
                                                                                    unit="kHz", scale=1, ndecimals=3
                                                                                ), group=self.name)

        # subsequences
        self.tickle_subsequence =                                               TickleDDS(self)

        # get relevant devices
        self.setattr_device('dds_modulation')

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # convert QLMS modulation to machine units
        self.freq_qlms_rabi_ftw_list =                                          np.array([
                                                                                    self.dds_modulation.frequency_to_ftw(freq_khz * kHz)
                                                                                    for freq_khz in self.freq_qlms_rabi_khz_list
                                                                                ])

        # create an array of values for the experiment to sweep
        # (i.e. DDS tickle frequency & readout FTW)
        self.config_qlms_rabi_list =                                            np.stack(np.meshgrid(self.freq_qlms_rabi_ftw_list, self.freq_readout_ftw_list), -1).reshape(-1, 2)
        np.random.shuffle(self.config_qlms_rabi_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_qlms_rabi_ftw_list) * len(self.freq_readout_ftw_list),
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.tickle_subsequence.record_dma()

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

            # sweep experiment config: heating time and readout frequency
            for config_vals in self.config_qlms_rabi_list:

                # extract values from config list
                freq_qlms_ftw =     config_vals[0]
                freq_readout_ftw =  config_vals[1]
                self.core.break_realtime()

                # set QLMS modulation frequency
                self.dds_modulation.set_mu(freq_qlms_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0)

                # set readout frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # QLMS tickle
                self.tickle_subsequence.run_dma()

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), freq_qlms_ftw)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)
