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
        self.setattr_argument("time_heating_rate_ms_list",                      PYONValue([1, 2, 5, 10, 50]))

        # run regular sideband cooling build
        super().build_experiment()

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # convert heating rate timings to machine units
        self.time_heating_rate_mu_list =                                        np.array([self.core.seconds_to_mu(time_ms * ms)
                                                                                          for time_ms in self.time_heating_rate_ms_list],
                                                                                         dtype=np.int64)

        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_heating_rate_list =                                         np.stack(np.meshgrid(self.time_heating_rate_mu_list, self.freq_readout_ftw_list), -1).reshape(-1, 2)
        self.config_heating_rate_list =                                         np.array(self.config_heating_rate_list, dtype=np.int64)
        np.random.shuffle(self.config_heating_rate_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_heating_rate_list),
                3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep experiment config: heating time and readout frequency
            for config_vals in self.config_heating_rate_list:

                # extract values from config list
                time_heating_delay_mu =     config_vals[0]
                freq_readout_ftw =          np.int32(config_vals[1])
                self.core.break_realtime()

                # set frequency
                self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # wait time to measure heating rate
                delay_mu(time_heating_delay_mu)

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(), time_heating_delay_mu)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)


    # ANALYSIS
    def analyze(self):
        """
        Fit resultant spectra with a sinc profile to extract n,
        then fit a line to extract heating rate
        """
        pass
