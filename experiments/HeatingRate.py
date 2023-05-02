import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class HeatingRate(SidebandCooling.SidebandCooling):
    """
    Experiment: Heating Rate

    Measures the heating rate by doing sideband cooling, then waiting
    a variable amount of time before readout.
    """
    name = 'Heating Rate'


    def build_experiment(self):
        # heating rate wait times
        self.setattr_argument("time_heating_rate_ms_list",                      PYONValue([1, 2]))

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list =                                        np.array([seconds_to_mu(time_ms * ms)
                                                                                          for time_ms in self.time_heating_rate_ms_list],
                                                                                         dtype=np.int64)

        # run preparations for sideband cooling
        super().prepare_experiment()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_heating_rate_mu_list) * len(self.freq_readout_ftw_list),
                3)


    # MAIN SEQUENCE
    @kernel
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

                    # wait time to measure heating rate
                    # delay_mu(time_heating_delay_mu)

                    # custom SBC readout
                    self.core_dma.playback_handle(_handle_sbc_readout)

                    # update dataset
                    with parallel:
                        self.update_results(freq_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
                        self.core.break_realtime()
