import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import TickleDDS
import LAX_exp.experiments.diagnostics.SidebandCooling as SidebandCooling


class QLMSRamsey(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Ramsey

    Quantum Logic Mass Spectroscopy - Ramsey Spectroscopy
    Cool the ions to the ground state of motion via sideband cooling,
    then conduct ramsey spectroscopy on a coherent state (created via tickle) to be read out via RSB/BSB comparison.
    """
    name = 'QLMSRamsey'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # QLMS configuration
        self.setattr_argument("time_qlms_ramsey_delay_ms_list",                Scannable(
                                                                                    default=[
                                                                                        RangeScan(1, 10, 11, randomize=True),
                                                                                        ExplicitScan([1])
                                                                                    ],
                                                                                    global_min=0.000001, global_max=1000000, global_step=1,
                                                                                    unit="ms", scale=1, ndecimals=3
                                                                                ), group=self.name)
        self.setattr_argument("freq_qlms_ramsey_khz_list",                      Scannable(
                                                                                    default=CenterScan(1101, 10, 0.5, randomize=True),
                                                                                    global_min=0, global_max=10000, global_step=1,
                                                                                    unit="kHz", scale=1, ndecimals=3
                                                                                ), group=self.name)

        # subsequences & devices
        self.tickle_subsequence =                                               TickleDDS(self)
        self.setattr_device('dds_parametric')

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()
        self.freq_sideband_readout_ftw_list =                                   self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list

        # convert QLMS parameters to machine units
        self.time_qlms_ramsey_delay_mu_list =                                   np.array([self.core.seconds_to_mu(time_ms * ms)
                                                                                          for time_ms in self.time_qlms_ramsey_delay_ms_list],
                                                                                         dtype=np.int64)
        self.freq_qlms_ramsey_ftw_list =                                        np.array([self.dds_parametric.frequency_to_ftw(freq_khz * kHz)
                                                                                          for freq_khz in self.freq_qlms_ramsey_khz_list])

        # create an array of values for the experiment to sweep
        # (i.e. DDS tickle frequency & readout FTW)
        self.config_qlms_ramsey_list =                                          np.stack(np.meshgrid(self.time_qlms_ramsey_delay_mu_list,
                                                                                                     self.freq_qlms_ramsey_ftw_list,
                                                                                                     self.freq_sideband_readout_ftw_list),
                                                                                         -1).reshape(-1, 3)
        self.config_qlms_ramsey_list =                                          np.array(self.config_qlms_ramsey_list, dtype=np.int64)
        np.random.shuffle(self.config_qlms_ramsey_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_qlms_ramsey_list),
                3)



    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):
            for config_vals in self.config_qlms_ramsey_list:

                # extract values from config list
                time_qlms_mu =      config_vals[0]
                freq_qlms_ftw =     np.int32(config_vals[1])
                freq_readout_ftw =  np.int32(config_vals[2])
                self.core.break_realtime()

                # prepare devices waveforms
                self.dds_parametric.set_mu(freq_qlms_ftw, asf=self.dds_parametric.ampl_modulation_asf, profile=0)
                self.qubit.set_mu(freq_readout_ftw, self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # QLMS ramsey
                self.tickle_subsequence.run()
                delay_mu(time_qlms_mu)
                self.tickle_subsequence.run()

                # sideband readout
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), freq_qlms_ftw, time_qlms_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()
