import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, ReadoutAdaptive, RescueIon


class CalibrationAdaptiveReadout(LAXExperiment, Experiment):
    """
    Calibration: Adaptive Readout

    Test adaptive, MLE-based readout (assuming a single-ion system).
    Technique from Alice H. Burrell thesis (2010, Lucas/Oxford).
    """
    name = 'Calibration Adaptive Readout'
    kernel_invariants = {
        # hardware objects & parameters
        'initialize_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'freq_qubit_ftw', 'att_qubit_mu', 'time_qubit_mu', 'profile_729_target'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=10000, precision=0, step=1, min=1, max=100000000))

        # sequence configuration
        self.setattr_argument("enable_doppler", BooleanValue(default=False))
        self.setattr_argument("enable_qubit",   BooleanValue(default=True))

        # adaptive readout parameters
        self.setattr_argument("time_bin_us",        NumberValue(default=8, precision=3, step=5, min=0.1, max=10000), group="readout")
        self.setattr_argument("error_threshold",    NumberValue(default=1e-2, precision=8, step=1e-2, min=1e-10, max=1.), group="readout")

        # qubit parameters
        self.setattr_argument("freq_qubit_mhz", NumberValue(default=101.1072, precision=6, step=1, min=50., max=400.), group="qubit")
        self.setattr_argument("time_qubit_us",  NumberValue(default=1.5, precision=2, step=5, min=0.1, max=10000), group="qubit")
        self.setattr_argument("att_qubit_db",   NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5), group="qubit")

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

        # subsequences
        self.profile_729_target =   6
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      ReadoutAdaptive(self, time_bin_us=self.time_bin_us,
                                                        error_threshold=self.error_threshold)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # qubit values
        self.freq_qubit_ftw =   self.qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)
        self.time_qubit_mu =    self.core.seconds_to_mu(self.time_qubit_us * us)
        self.att_qubit_mu =     att_to_mu(self.att_qubit_db * dB)

        # '''
        # CREATE EXPERIMENT CONFIG
        # '''
        # # create an array of values for the experiment to sweep
        # self.config_experiment_list = np.stack(np.meshgrid(
        #     freq_deshelve_ftw_list,
        #     ampl_deshelve_pct_list,
        #     time_deshelve_mu_list
        # ), -1).reshape(-1, 3)
        # self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        # np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions, 3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.core.break_realtime()

        # configure qubit beam
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit.set_mu(self.freq_qubit_ftw, asf=self.qubit.ampl_qubit_asf,
                          profile=self.profile_729_target)
        self.qubit.set_profile(self.profile_729_target)
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # main sequence
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # initialize ion
            if self.enable_doppler or self.enable_qubit:
                self.initialize_subsequence.run_dma()

            # apply qubit pulse
            if self.enable_qubit:
                self.qubit.on()
                delay_mu(self.time_qubit_mu)
                self.qubit.off()

            # adaptive readout
            results = self.readout_subsequence.run()

            # finish up and add slack
            self.update_results(results[0], results[1], results[2])
            delay_mu(50000)

            # periodically check termination
            if trial_num % 100 == 1:
                if self.scheduler.check_termination():
                    self.core.break_realtime()
                    break

            # resuscitate ion
            self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

    def analyze_experiment(self):
        """
        Print summary statistics.
        """
        print("############# RESULT SUMMARY #############")
        # collate data and ensure correct shape for processing
        data_bright =   self.results[self.results[:, 0] == 1]
        data_dark =     self.results[self.results[:, 0] == 0]
        data_idk =      self.results[self.results[:, 0] == -1]
        if len(data_bright) == 0:   data_bright = np.ones((1, np.shape(self.results)[1])) * np.nan
        if len(data_dark) == 0:     data_dark = np.ones((1, np.shape(self.results)[1])) * np.nan
        if len(data_idk) == 0:      data_idk = np.ones((1, np.shape(self.results)[1])) * np.nan

        # calculate summary statistics and store as dataset
        num_bright, num_dark, num_idk = (len(data_bright[:, 0]), len(data_dark[:, 0]), len(data_idk[:, 0]))
        det_time_bright, det_time_std_bright, det_time_dark, det_time_std_dark = (
            self.core.mu_to_seconds(np.mean(data_bright[:, 2])), self.core.mu_to_seconds(np.std(data_bright[:, 2])),
            self.core.mu_to_seconds(np.mean(data_dark[:, 2])), self.core.mu_to_seconds(np.std(data_dark[:, 2]))
        )
        rates_bright, rates_dark, rates_idk = tuple(
            arr[:, 1] / (self.core.mu_to_seconds(arr[:, 2])) * 3e-3
            for arr in (data_bright, data_dark, data_idk)
        )
        self.set_dataset("det_time_s_bright_dark", [det_time_bright, det_time_dark])
        self.set_dataset("det_time_std_s_bright_dark", [det_time_std_bright, det_time_std_dark])
        self.set_dataset("count_rate_3ms_bright_dark_idk", [np.mean(rates_bright), np.mean(rates_dark), np.mean(rates_idk)])
        self.set_dataset("count_rate_std_3ms_bright_dark_idk", [np.std(rates_bright), np.std(rates_dark), np.std(rates_idk)])

        print("Discrimination Results (%, total events):"
              "\n\tBright:\t\t{:.3f}% ({:d})\n\tDark:\t\t{:.3f}% ({:d})\n\tIndeterminate:\t{:.3f}% ({:d})\n".format(
            num_bright / self.repetitions * 100., num_bright,
            num_dark / self.repetitions * 100., num_dark,
            num_idk / self.repetitions * 100., num_idk
        ))
        print("Time to Detection (us):"
              "\n\tBright:\t\t{:.1f} +/- {:.3g}\n\tDark:\t\t{:.1f}  +/- {:.3g}\n".format(
            det_time_bright / us, det_time_std_bright / us, det_time_dark / us, det_time_std_dark / us
        ))
        print("Count Rates (per 3ms):"
              "\n\tBright:\t\t{:.2f} +/- {:.3g}\n\tDark:\t\t{:.2f} +/- {:.3g}\n\tIndeterminate:\t{:.2f} +/- {:.3g}\n".format(
            np.mean(rates_bright), np.std(rates_bright),
            np.mean(rates_dark), np.std(rates_dark),
            np.mean(rates_idk), np.std(rates_idk)
        ))
