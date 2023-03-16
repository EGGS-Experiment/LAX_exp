import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
import LAX_exp.experiments.SidebandCooling as SidebandCooling


class EGGSHeating(SidebandCooling.SidebandCooling):
    """
    Experiment: EGGS Heating

    Apply the EGGS heating RF immediately after sideband cooling
    and read out the final phonon number.
    """
    name = 'EGGS Heating'

    def build_experiment(self):
        # heating rate wait times
        self.setattr_argument("time_heating_rate_ms_list",                      PYONValue([1, 2]))

        # eggs heating
        self.setattr_argument("time_eggs_heating_ms",                       NumberValue(default=2, ndecimals=5, step=1, min=0.000001, max=10000))
        self.setattr_argument("att_eggs_heating_db",                        NumberValue(default=10, ndecimals=1, step=0.5, min=10, max=31.5))
        self.setattr_argument("freq_eggs_heating_mhz_carrier_list",         Scannable(
                                                                                default=CenterScan(85.1, 0.020, 0.001, randomize=True),
                                                                                global_min=30, global_max=400, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ))
        self.setattr_argument("freq_eggs_heating_secular_mhz_list",         Scannable(
                                                                                default=CenterScan(1.41, 0.001, 0.0001, randomize=True),
                                                                                global_min=0, global_max=10000, global_step=0.1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ))

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list =                                        np.array([seconds_to_mu(time_ms * ms)
                                                                                          for time_ms in self.time_heating_rate_ms_list], dtype=np.int64)

        # run regular sideband cooling prepare
        super().prepare_experiment()

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_heating_rate_mu_list) * len(self.freq_readout_ftw_list),
                3)


    # MAIN SEQUENCE
    @kernel
    def run_main(self):
        self.core.reset()

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

                    # wait given time
                    delay_mu(time_heating_delay_mu)

                    # set readout profile for qubit
                    self.qubit.carrier()

                    # rabi flop
                    self.rabiflop_subsequence.run_dma()

                    # read out
                    self.readout_subsequence.run_dma()

                    # update dataset
                    with parallel:
                        self.update_results(freq_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
                        self.core.break_realtime()
