from artiq.experiment import *


class ParameterSet(EnvExperiment):
    """
    Parameter Set
    Sets all relevant parameters for experiment.
    Important parameters are supposed to persist; this is for problems/resets.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """

        # devices
        device_parameters = {
            # PMT
            "pmt_input_channel":            1,
            "pmt_gating_edge":              "rising",

            # DDS - boards & channels
            "dds_board_num":                1,
            "dds_probe_channel":            0,
            "dds_pump_channel":             1,
            "dds_repump_cooling_channel":   2,
            "dds_repump_qubit_channel":     3,

            "dds_board_qubit_num":          0,
            "dds_qubit_channel":            0,
        }

        # timing
        timing_parameters = {
            # preparation
            "time_repump_qubit_us":         100,
            "time_redist_us":               500,

            # cooling
            "time_doppler_cooling_us":      100,

            # readout
            "time_readout_us":              100
        }

        # DDS waveforms
        dds_parameters = {
            # frequency
            "freq_probe_mhz":           90,
            "freq_pump_cooling_mhz":    90,
            "freq_pump_readout_mhz":    92,
            "freq_repump_cooling_mhz":  110,
            "freq_repump_qubit_mhz":    110,
            "freq_qubit_mhz":           110.771,

            # amplitude
            "ampl_probe_pct":           50,
            "ampl_pump_pct":            50,
            "ampl_repump_cooling_pct":  50,
            "ampl_repump_qubit_pct":    50,
            "ampl_qubit_pct":           50,

            # attenuation
            "att_probe_dB":             21,
            "att_pump_dB":              22,
            "att_repump_cooling_dB":    22,
            "att_repump_qubit_dB":      22,
            "att_qubit_dB":             18
        }

        # consolidate parameter dictionaries
        consolidated_parameters = {
            **device_parameters,
            **timing_parameters,
            **dds_parameters
        }

        # set parameters as global persistent datasets
        for parameter_name, parameter_value in consolidated_parameters.items():
            self.set_dataset(parameter_name, parameter_value, broadcast=True, persist=True)


    def run(self):
        pass

