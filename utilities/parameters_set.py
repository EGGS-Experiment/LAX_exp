from artiq.experiment import *


class ParameterSet(EnvExperiment):
    """
    Parameter Set
    Sets all relevant parameters (persistent) for experiment.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """

        # readout parameters
        readout_parameters = {
            # PMT
            "pmt.pmt_input_channel":                0,
            "pmt.pmt_gating_edge":                  "rising",
            "pmt.pmt_discrimination":               11,

            # photodiode
            "photodiode.photodiode_channel":        0,
            "photodiode.photodiode_gain":           1
        }


        # rf parameters
        rf_parameters = {
            "rf.ttl_channel_rf_modulation":         3
        }


        # beam parameters
        beam_parameters = {
            # urukul boards
            "beams.dds_board_num":                  1,
            "beams.dds_board_qubit_num":            0,
            "beams.dds_board_tickle_num":           0,

            # urukul channels
            "beams.dds_qubit_channel":              1,
            "beams.dds_probe_channel":              0,
            "beams.dds_pump_channel":               1,
            "beams.dds_repump_cooling_channel":     2,
            "beams.dds_repump_qubit_channel":       3,

            # frequency
            "beams.freq_probe_redist_mhz":          110,
            "beams.freq_pump_cooling_mhz":          110,
            "beams.freq_pump_readout_mhz":          110,
            "beams.freq_repump_cooling_mhz":        110,
            "beams.freq_repump_qubit_mhz":          110,
            "beams.freq_qubit_mhz":                 110.771,

            # amplitude
            "beams.ampl_probe_redist_pct":          50,
            "beams.ampl_pump_pct":                  50,
            "beams.ampl_pump_cooling_pct":          50,
            "beams.ampl_pump_readout_pct":          50,
            "beams.ampl_repump_cooling_pct":        50,
            "beams.ampl_repump_qubit_pct":          50,
            "beams.ampl_qubit_pct":                 50
        }


        # external devices
        external_parameters = {
            # devices
            "function_generator.ttl_channel":       9,
        }


        # timing
        timing_parameters = {
            # preparation
            "time_repump_qubit_us":                 100,
            "time_redist_us":                       500,
            "time_probe_us":                        50,

            # cooling
            "time_doppler_cooling_us":              800,

            # readout
            "time_readout_us":                      500,

            # general
            "time_profileswitch_delay_us":          1
        }


        # consolidate parameter dictionaries
        consolidated_parameters = {
            **readout_parameters,
            **rf_parameters,
            **beam_parameters,
            **external_parameters,
            **timing_parameters
        }

        # set parameters as global persistent datasets
        for parameter_name, parameter_value in consolidated_parameters.items():
            self.set_dataset(parameter_name, parameter_value, broadcast=True, persist=True)


    def run(self):
        pass
