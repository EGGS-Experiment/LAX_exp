import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
import LAX_exp.experiments.diagnostics.LaserScan as LaserScan


class LaserScanMulti(LaserScan.LaserScan):
    """
    Experiment: Laser Scan Multi

    Does a 729nm laser scan; resets the ion(s) every shot.
    Supports sine-squared pulse shaping.
    Scans multiple user-selectable ranges.
    """
    name = 'Laser Scan Multi'
    # note: no need to specify kernel invariants since parent specifies all of them for us

    def build_experiment(self):
        # call parent build
        super().build_experiment()

        # laser scan multi - arguments
        self.setattr_argument("scan_type", EnumerationValue(["Scan", "Scan+Sideband1", "Scan+Sideband2", "Scan+Sideband3",
                                                             "Scan+1+2", "Scan+All", "Sideband1", "Sideband2", "Sideband3",
                                                             "1+2", "All Sidebands"],
                                                            default="Scan+All"), group="Multiscan")
        self.setattr_argument("freq_sideband_1_khz",    NumberValue(default=702.1, min=-4e5, max=4e5, step=1, unit="kHz", scale=1, precision=3), group="Multiscan")
        self.setattr_argument("freq_sideband_2_khz",    NumberValue(default=1303.1, min=-4e5, max=4e5, step=10, unit="kHz", scale=1, precision=3), group="Multiscan")
        self.setattr_argument("freq_sideband_3_khz",    NumberValue(default=1592.1, min=-4e5, max=4e5, step=1, unit="kHz", scale=1, precision=3), group="Multiscan")

    def prepare_experiment(self):
        # pre-convert variables for convenience
        carrier_array_mhz = np.array(list(self.freq_qubit_scan_mhz))
        freq_sideband_arrays_aom_mhz_list = np.array([
            carrier_array_mhz - self.freq_sideband_1_khz * (kHz/MHz) / 2.,
            carrier_array_mhz - self.freq_sideband_2_khz * (kHz/MHz) / 2.,
            carrier_array_mhz - self.freq_sideband_3_khz * (kHz/MHz) / 2.
        ])

        # process scan type: only need to update self.freq_qubit_scan_mhz since parent does
        # all processing based on that variable
        if self.scan_type == "Scan":
            self.freq_qubit_scan_mhz = carrier_array_mhz
        elif self.scan_type == "Scan+Sideband1":
            self.freq_qubit_scan_mhz = np.concatenate((carrier_array_mhz,
                                                      freq_sideband_arrays_aom_mhz_list[0]))
        elif self.scan_type == "Scan+Sideband2":
            self.freq_qubit_scan_mhz = np.concatenate((carrier_array_mhz,
                                                      freq_sideband_arrays_aom_mhz_list[1]))
        elif self.scan_type == "Scan+Sideband3":
            self.freq_qubit_scan_mhz = np.concatenate((carrier_array_mhz,
                                                      freq_sideband_arrays_aom_mhz_list[2]))
        elif self.scan_type == "Scan+1+2":
            self.freq_qubit_scan_mhz = np.concatenate((carrier_array_mhz,
                                                      freq_sideband_arrays_aom_mhz_list[0],
                                                      freq_sideband_arrays_aom_mhz_list[1]))
        elif self.scan_type == "Scan+All":
            self.freq_qubit_scan_mhz = np.concatenate((carrier_array_mhz,
                                                      freq_sideband_arrays_aom_mhz_list[0],
                                                      freq_sideband_arrays_aom_mhz_list[1],
                                                      freq_sideband_arrays_aom_mhz_list[2]))
        elif self.scan_type == "Sideband1":
            self.freq_qubit_scan_mhz = freq_sideband_arrays_aom_mhz_list[0]
        elif self.scan_type == "Sideband2":
            self.freq_qubit_scan_mhz = freq_sideband_arrays_aom_mhz_list[1]
        elif self.scan_type == "Sideband3":
            self.freq_qubit_scan_mhz = freq_sideband_arrays_aom_mhz_list[2]
        elif self.scan_type == "1+2":
            self.freq_qubit_scan_mhz = np.concatenate((freq_sideband_arrays_aom_mhz_list[0],
                                                       freq_sideband_arrays_aom_mhz_list[1]))
        elif self.scan_type == "All":
            self.freq_qubit_scan_mhz = np.concatenate((freq_sideband_arrays_aom_mhz_list[0],
                                                       freq_sideband_arrays_aom_mhz_list[1],
                                                       freq_sideband_arrays_aom_mhz_list[2]))

        # call parent prepare_experiment
        super().prepare_experiment()

