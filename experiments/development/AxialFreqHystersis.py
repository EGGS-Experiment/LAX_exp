import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.LaserScan as LaserScan


class AxialFreqHystersis(LaserScan.LaserScan):
    """
    Experiment: Axial Freq Hystersis

    Does a 729nm laser scan; resets the ion(s) every shot.
    Supports sine-squared pulse shaping.
    Scans DC output for axial trapping potential to determine hystersis
    """
    name = 'Axial Freq Hystersis'
    # note: no need to specify kernel invariants since parent specifies all of them for us

    def build_experiment(self):
        # call parent build
        super().build_experiment()

        self.setattr_argument('dc_voltage_endcap_scan', Scannable([
            ExplicitScan([-1,0,1]),
            RangeScan(-1,1, 10, randomize=False),
            CenterScan(0,5,1, randomize=False)
        ],
            global_max=400,
            global_min=0, global_step=0.01,
            unit = 'V', precision=2
        ), group='voltage_scan',
        tooltip='relative change in voltage to BOTH endcaps'
                'NOTE BENE: these values are how much the endcap voltages are CHANGED by, they are'
                'NOT THE ABSOLUTE voltages')

        self.setattr_argument("num_voltage_sweeps", NumberValue(1,
                                                                step=1, precision=0,
                                                                max = 100, min=1),
                              group='voltage_scan',
                              tooltip='Number of voltage sweeps to perform when determining if there is hystersis')

        # laser scan multi - arguments
        self.setattr_argument('trap_dc')

    def prepare_experiment(self):
        super().prepare_experiment()

        west_endcap_voltage = self.trap_dc.get_west_endcap_voltage()
        east_endcap_voltage = self.trap_dc.get_west_endcap_voltage()

        self.east_endcap_voltage_scan = east_endcap_voltage + array(self.dc_voltage_endcap_scan)
        self.west_endcap_voltage_scan = west_endcap_voltage + array(self.dc_voltage_endcap_scan)

        self.endcap_voltage_scan = zip(self.east_endcap_voltage_scan, self.west_endcap_voltage_scan,
                                       array(self.dc_voltage_scan)
                                       )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)

    def run_main(self) -> TNone:
        for num in range(self.num_voltage_sweeps):
            for voltages in self.endcap_voltage_scan:

                self.trap_dc.set_east_endcap_voltage(voltages[0])
                self.trap_dc.get_west_endcap_voltage(voltages[1])

                delay_mu(1e9)

                for trial_num in range(self.repetitions):
                    for config_vals in self.config_experiment_list:

                        ### PREPARE & CONFIGURE ###
                        # extract values from config list
                        freq_ftw = int32(config_vals[0])
                        time_holdoff_mu = config_vals[1]

                        # prepare relevant beams
                        self.core.break_realtime()
                        if self.enable_pulseshaping:
                            self.qubit.set_ftw(freq_ftw)
                        else:
                            self.qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf,
                                              profile=self.profile_729_readout,
                                              phase_mode=PHASE_MODE_CONTINUOUS)
                        delay_mu(10000)

                        # wait for linetrigger
                        if self.enable_linetrigger:
                            self.trigger_line.trigger(self.trigger_line.time_timeout_mu, time_holdoff_mu)

                        ### MAIN SHOT ###
                        # initialize ion in S-1/2 state
                        self.initialize_subsequence.run_dma()

                        # fire spectroscopy pulse
                        if self.enable_pulseshaping:
                            self.pulseshape_subsequence.run()
                        else:
                            self.qubit.on()
                            delay_mu(self.time_qubit_mu)
                            self.qubit.off()

                        # read out counts & clean up loop
                        self.readout_subsequence.run_dma()
                        self.rescue_subsequence.resuscitate()
                        self.initialize_subsequence.slack_rescue()
                        counts = self.readout_subsequence.fetch_count()

                        # store results in dataset
                        self.rescue_subsequence.detect_death(counts)
                        self.update_results(freq_ftw, counts, time_holdoff_mu, voltages[2])

                        # rescue ion as needed & support graceful termination
                        self.core.break_realtime()
                        self.rescue_subsequence.run(trial_num)
                        self.check_termination()

    def analyze_experiment(self):
        pass