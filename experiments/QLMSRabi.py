import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class QLMSRabi(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Rabi

    Quantum Logic Mass Spectroscopy - Rabi Flopping
    Cool the ions to the ground state of motion via sideband cooling,
    then apply a quadrupole tickle to create a coherent state to be read out via RSB/BSB comparison.
    """
    name = 'QLMSRabi'


    def build_experiment(self):
        # QLMS configuration
        self.setattr_argument("time_qlms_heating_ms",                           NumberValue(default=2, ndecimals=5, step=1, min=0.000001, max=10000), group=self.name)
        self.setattr_argument("att_qlms_heating_db",                            NumberValue(default=30, ndecimals=1, step=0.5, min=0, max=31.5), group=self.name)
        self.setattr_argument("freq_qlms_heating_khz_list",                     Scannable(
                                                                                    default=CenterScan(1558, 100, 1, randomize=True),
                                                                                    global_min=0, global_max=10000, global_step=1,
                                                                                    unit="MHz", scale=1, ndecimals=3
                                                                                ))
        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # todo: convert QLMS heating to machine units

        # run preparations for sideband cooling
        super().prepare_experiment()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_heating_rate_mu_list) * len(self.freq_readout_ftw_list),
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep times to measure heating rate
            for time_heating_delay_mu in self.time_heating_rate_mu_list:

                # sweep frequency
                for freq_ftw in self.freq_readout_ftw_list:

                    # set frequency
                    self.qubit.set_mu(freq_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                    self.core.break_realtime()

                    # initialize ion in S-1/2 state
                    self.initialize_subsequence.run_dma()

                    # sideband cool
                    self.sidebandcool_subsequence.run_dma()

                    # tmp remove
                    self.core.break_realtime()
                    self.dds_mod.set(sideband_freq_hz, amplitude=0.35, profile=0)
                    self.dds_mod.cpld.set_profile(0)
                    self.dds_mod.cpld.io_update.pulse_mu(8)
                    self.dds_mod.set_att(self.att_eggs_heating_db)
                    self.dds_mod.cfg_sw(True)
                    delay_mu(self.time_eggs_heating_mu)
                    self.dds_mod.cfg_sw(False)
                    # tmp remove

                    # custom SBC readout
                    self.core_dma.playback_handle(_handle_sbc_readout)

                    # update dataset
                    with parallel:
                        self.update_results(freq_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
                        self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)
