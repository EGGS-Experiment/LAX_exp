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
            "pmt.pmt_input_channel":                                    0,
            "pmt.pmt_gating_edge":                                      "rising",
            "pmt.pmt_discrimination":                                   11,

            # photodiode
            "photodiode.photodiode_channel":                            0,
            "photodiode.photodiode_gain":                               1
        }


        # rf parameters
        rf_parameters = {
            "rf.ttl_channel_rf_modulation":                             3
        }


        # beam parameters
        beam_parameters = {
            # urukul boards
            "beams.dds_board.dds_board_num":                            1,
            "beams.dds_board.dds_board_qubit_num":                      0,
            "beams.dds_board.dds_board_tickle_num":                     0,

            # urukul channels
            "beams.dds_channel.dds_qubit_channel":                      1,
            "beams.dds_channel.dds_probe_channel":                      0,
            "beams.dds_channel.dds_pump_channel":                       1,
            "beams.dds_channel.dds_repump_cooling_channel":             2,
            "beams.dds_channel.dds_repump_qubit_channel":               3,

            # frequency
            "beams.freq_mhz.freq_probe_redist_mhz":                     110.0,
            "beams.freq_mhz.freq_probe_spinpol_mhz":                    110.0,
            "beams.freq_mhz.freq_pump_cooling_mhz":                     106.0,
            "beams.freq_mhz.freq_pump_readout_mhz":                     110.0,
            "beams.freq_mhz.freq_repump_cooling_mhz":                   110.0,
            "beams.freq_mhz.freq_repump_qubit_mhz":                     110.0,
            "beams.freq_mhz.freq_qubit_carrier_mhz":                    110.771,
            "beams.freq_mhz.freq_qubit_rsb_mhz":                        104.012,
            "beams.freq_mhz.freq_qubit_bsb_mhz":                        105.214,

            # amplitude
            "beams.ampl_pct.ampl_probe_redist_pct":                     16.0,
            "beams.ampl_pct.ampl_probe_spinpol_pct":                    50.0,
            "beams.ampl_pct.ampl_pump_cooling_pct":                     21.0,
            "beams.ampl_pct.ampl_pump_readout_pct":                     45.0,
            "beams.ampl_pct.ampl_repump_cooling_pct":                   17.0,
            "beams.ampl_pct.ampl_repump_qubit_pct":                     15.0,
            "beams.ampl_pct.ampl_qubit_pct":                            50.0,
            'beams.ampl_pct.ampl_tickle_radial_pct':                    50.0,
            'beams.ampl_pct.ampl_tickle_radial_pct':                    50.0
        }


        # motion
        motion_parameters = {
            # secular frequencies
            "motion.freq_khz.freq_secular_x_khz":                       1450.0,
            "motion.freq_khz.freq_secular_y_khz":                       1375.0,
            "motion.freq_khz.freq_secular_z_khz":                       550.0,

            # amplitude
            "beams.ampl_pct.ampl_probe_redist_pct":                     50.0,
            "beams.ampl_pct.ampl_probe_spinpol_pct":                    50.0,
            "beams.ampl_pct.ampl_pump_pct":                             50.0,
            "beams.ampl_pct.ampl_pump_cooling_pct":                     50.0,
            "beams.ampl_pct.ampl_pump_readout_pct":                     50.0,
            "beams.ampl_pct.ampl_repump_cooling_pct":                   21.0,
            "beams.ampl_pct.ampl_repump_qubit_pct":                     50.0,
            "beams.ampl_pct.ampl_qubit_pct":                            50.0
        }

        # external devices
        external_parameters = {
            # devices
            "function_generator.ttl_channel_function_generator":        9,
        }


        # timing
        timing_parameters = {
            # state preparation
            "timing.time_repump_qubit_us":                              100,
            "timing.time_redist_us":                                    500,
            "timing.time_spinpol_us":                                   500,
            "timing.time_probe_us":                                     50,

            # cooling & readout
            "timing.time_doppler_cooling_us":                           800,
            "timing.time_readout_us":                                   500,
            "timing.time_rabiflop_us":                                  50,

            # tickle
            "timing.time_tickle_us":                                    500,

            # general
            "timing.time_profileswitch_delay_us":                       1
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
