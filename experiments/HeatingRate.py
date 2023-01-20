import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.experiments.SidebandCooling import SidebandCooling2


class HeatingRate(SidebandCooling2, Experiment):
    """
    Experiment: Heating Rate

    Measures the heating rate by doing sideband cooling, then waiting
    a variable amount of time before readout.
    """

    name = 'Heating Rate'

    def build_experiment(self):
        # heating rate wait times
        self.setattr_argument("time_heating_rate_ms_list",              PYONValue([1, 2]))

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list =                                np.array(seconds_to_mu(self.time_heating_rate_mu_list * ms), dtype=np.int64)

        # run regular sideband cooling prepare
        super().prepare_experiment()
        # todo: think about how to do dataset setup since prepare sets up datasets

    # MAIN SEQUENCE
    @kernel
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # sweep times to measure heating rate
            for time_heating_delay_mu in self.time_heating_rate_list_mu:

                # sweep frequency
                for freq_ftw in self.freq_qubit_scan_ftw:

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
                        self.update_dataset(freq_ftw, self.readout_subsequence.fetch_count())
                        self.core.break_realtime()

            self.set_dataset('management.completion_pct', (trial_num + 1) / self.repetitions * 100., broadcast=True, persist=True, archive=False)
