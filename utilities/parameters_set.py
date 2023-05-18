from artiq.experiment import *


class ParameterSet(EnvExperiment):
    """
    Parameter Set
    Sets all relevant parameters (persistent) for experiment.
    """

    def prepare(self):
        """
        Set devices and arguments for the experiment.
        """

        # readout parameters
        readout_parameters = {
            # PMT
            "pmt.gating_edge":                                          "rising",
            "pmt.time_pmt_gating_us":                                   100,

            # todo: distinguish between the 3 photodiodes
            # photodiode
            "photodiode.sampler_channel":                               0,
            "photodiode.gain":                                          1
        }

        # todo: maybe merge pmt, linetrigger, and rf parameters into one holding block?
        # ttl parameters
        ttl_parameters = {
            # linetrigger
            "linetrigger.gating_edge":                                  'rising',
            "linetrigger.time_timeout_ms":                              100.0,
            "linetrigger.time_holdoff_us":                              100.0
        }


        # rf parameters
        rf_parameters = {
            "rf.gating_edge":                                           "rising",
            "rf.time_rf_holdoff_us":                                    100.0,
            "rf.time_rf_gating_ns":                                     150.0,

        }


        # dds parameters
        dds_parameters = {
            # amplitude
            "dds.ampl_pct.ampl_modulation_pct":                         35.0,
        }


        # beam parameters
        beam_parameters = {
            # frequency
            "beams.freq_mhz.freq_probe_spinpol_mhz":                    114.0,
            "beams.freq_mhz.freq_pump_cooling_mhz":                     110.0,
            "beams.freq_mhz.freq_pump_readout_mhz":                     110.0,
            "beams.freq_mhz.freq_pump_rescue_mhz":                      110.0,
            "beams.freq_mhz.freq_repump_cooling_mhz":                   110.0,
            "beams.freq_mhz.freq_repump_qubit_mhz":                     110.0,
            "beams.freq_mhz.freq_qubit_carrier_mhz":                    110.771,
            "beams.freq_mhz.freq_qubit_rsb_mhz":                        104.012,
            "beams.freq_mhz.freq_qubit_bsb_mhz":                        105.214,

            # amplitude
            "beams.ampl_pct.ampl_probe_spinpol_pct":                    20.0,
            "beams.ampl_pct.ampl_pump_cooling_pct":                     19.0,
            "beams.ampl_pct.ampl_pump_readout_pct":                     45.0,
            "beams.ampl_pct.ampl_pump_rescue_pct":                      45.0,
            "beams.ampl_pct.ampl_repump_cooling_pct":                   15.0,
            "beams.ampl_pct.ampl_repump_qubit_pct":                     15.0,
            "beams.ampl_pct.ampl_qubit_pct":                            50.0,
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
            "timing.time_rescue_us":                                    1,

            # standard
            "timing.time_spinpol_us":                                   20,
            "timing.time_doppler_cooling_us":                           3000,
            "timing.time_readout_us":                                   500,
            "timing.time_probe_us":                                     50,
            "timing.time_repump_qubit_us":                              100,

            # qubit manipulation
            "timing.time_rabiflop_us":                                  50,
        }


        # calibration parameters
        calibration_parameters = {
            # qubit carrier
            "calibration.qubit.freq_carrier_mhz":                       104.364,
            "calibration.qubit.cavity_drift_rate_hz_per_s":             0.264,
            "calibration.qubit.calibration_timestamp":                  1673141207.260716,

            # secular frequency
            "calibration.secular.freq_secular_x_mhz":                   1.600,
            "calibration.secular.freq_secular_y_mhz":                   1.480,
            "calibration.secular.freq_secular_z_mhz":                   0.692,
            "calibration.secular.calibration_timestamp":                1673141207.260716,

            # pmt counts
            "calibration.pmt.sample_num":                               100,
            "calibration.pmt.counts_signal":                            17,
            "calibration.pmt.counts_noise":                             1.4,
            "calibration.pmt.calibration_timestamp":                    1673141207.260716
        }


        # experiment management parameters
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
            **ttl_parameters,
            **beam_parameters,
            **external_parameters,
            **timing_parameters,
            **management_parameters,
            **calibration_parameters
        }

        # set parameters as global persistent datasets
        for parameter_name, parameter_value in consolidated_parameters.items():
            self.set_dataset(parameter_name, parameter_value, broadcast=True, persist=True)


    def run(self):
        pass
