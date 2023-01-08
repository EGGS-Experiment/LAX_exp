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
            "photodiode.photodiode_channel_sampler":                    0,
            "photodiode.photodiode_gain":                               1
        }


        # rf parameters
        rf_parameters = {
            "rf.ttl_channel_rf_modulation":                             3
        }


        # beam parameters
        beam_parameters = {
            # frequency
            "beams.freq_mhz.freq_probe_spinpol_mhz":                    114.0,
            "beams.freq_mhz.freq_pump_cooling_mhz":                     110.0,
            "beams.freq_mhz.freq_pump_readout_mhz":                     110.0,
            "beams.freq_mhz.freq_repump_cooling_mhz":                   110.0,
            "beams.freq_mhz.freq_repump_qubit_mhz":                     110.0,
            "beams.freq_mhz.freq_qubit_carrier_mhz":                    110.771,
            "beams.freq_mhz.freq_qubit_rsb_mhz":                        104.012,
            "beams.freq_mhz.freq_qubit_bsb_mhz":                        105.214,

            # amplitude
            "beams.ampl_pct.ampl_probe_spinpol_pct":                    20.0,
            "beams.ampl_pct.ampl_pump_cooling_pct":                     19.0,
            "beams.ampl_pct.ampl_pump_readout_pct":                     45.0,
            "beams.ampl_pct.ampl_repump_cooling_pct":                   15.0,
            "beams.ampl_pct.ampl_repump_qubit_pct":                     15.0,
            "beams.ampl_pct.ampl_qubit_pct":                            50.0,
            'beams.ampl_pct.ampl_tickle_radial_pct':                    50.0
        }


        # external devices
        external_parameters = {
            # devices
            "function_generator.ttl_channel_function_generator":        9,
        }


        # timing
        timing_parameters = {
            # general
            "timing.time_profileswitch_delay_us":                       1,

            # standard
            "timing.time_spinpol_us":                                   20,
            "timing.time_doppler_cooling_us":                           3000,
            "timing.time_readout_us":                                   500,
            "timing.time_probe_us":                                     50,
            "timing.time_repump_qubit_us":                              100,

            # qubit manipulation
            "timing.time_rabiflop_us":                                  50,

            # tickle
            "timing.time_tickle_us":                                    500
        }


        # calibration parameters
        calibration_parameters = {
            # qubit carrier
            "calibration.qubit.freq_carrier_mhz":                       104.364,
            "calibration.qubit.cavity_drift_rate_hz_per_s":             0.264,
            "calibration.qubit.time_last_calibrated_timestamp":         1673141207.260716,

            # secular frequency
            "calibration.secular.freq_secular_x_mhz":                   1.600,
            "calibration.secular.freq_secular_x_mhz":                   1.480,
            "calibration.secular.freq_secular_x_mhz":                   0.692,

            # pmt counts
            "calibration.pmt.counts_signal":                            17,
            "calibration.pmt.counts_noise":                             1.4,
        }


        management_parameters = {
            # experiments
            "management.completion_pct":                                0.00,

            # datasets
            "management.dataset_save_locations":                        ["Z:\Motion\Data"]
        }

        # consolidate parameter dictionaries
        consolidated_parameters = {
            **readout_parameters,
            **rf_parameters,
            **beam_parameters,
            **external_parameters,
            **timing_parameters,
            **management_parameters
        }

        # set parameters as global persistent datasets
        for parameter_name, parameter_value in consolidated_parameters.items():
            self.set_dataset(parameter_name, parameter_value, broadcast=True, persist=True)


    def run(self):
        pass
