import numpy as np
from artiq.experiment import *
from LAX_exp.language import *

from copy import deepcopy
from random import shuffle

def update_deep(dict_dj, keyspec, val):
    """
    Modify values within a nested dict using keys in dot form.
    """
    keylist = keyspec.split('.')
    for key in keylist[:-1]:
        dict_dj = dict_dj[key]
    dict_dj[keylist[-1]] = val


class BatchSubmission(EnvExperiment):
    """
    idk document i guess
    """

    def build(self):
        self.scheduler = self.get_device("scheduler")

    def get_expid(self):
        expid = {
            "devarg_override": {},
            "log_level": 30,
            "file": "experiments\\ion_spectrum_analyzer\\SuperDuperResolutionAmpl.py",
            "class_name": "SuperDuperResolutionAmpl",
            "repo_rev": "34150853bfa4366dc0a345ae93b40f0f7e3b60e8",

            "arguments": {
                "repetitions": 200,
                "randomize_config": True,
                "sub_repetitions": 1,
                "readout_type": "RAP",


                "freq_rsb_readout_mhz_list": {
                    "center": 100.4596, "span": 0.03, "step": 0.0006,
                    "randomize": True, "seed": None, "ty": "CenterScan"
                },
                "freq_bsb_readout_mhz_list": {
                    "center": 101.7189, "span": 0.03, "step": 0.0006,
                    "randomize": True, "seed": None, "ty": "CenterScan"
                },
                "ampl_sideband_readout_pct": 50.0,
                "att_sideband_readout_db": 8.0,
                "time_sideband_readout_us": 46.3,
                "time_readout_us_list": {"sequence": [46.3], "ty": "ExplicitScan"},


                "enable_rescue": False,
                "enable_resuscitate": False,
                "enable_aperture": False,
                "add_397nm_spinpol": False,
                "repetitions_per_rescue": 1,


                "ampl_fock_pct": 50.0,
                "fock_num_prep": 1,
                "fock_num_read": "1",
                "enable_final_rap": True,


                "freq_phaser_carrier_mhz_list": {
                    "sequence": [86.0], "ty": "ExplicitScan"
                },
                "freq_sweep_arr": "[-1.0, 1.0, 0.0, 0.0]",
                "freq_osc_sweep_khz_list": {
                    "center": 0.0, "span": 20.0, "step": 0.5,
                    "randomize": True, "seed": None, "ty": "CenterScan"
                },

                "phase_sweep_arr": "[1.0, 0.0, 0.0, 0.0]",
                "phase_osc_sweep_turns_list": {"sequence": [0.0], "ty": "ExplicitScan"},
                "phase_global_ch1_turns_list": {"sequence": [0.1487], "ty": "ExplicitScan"},


                "enable_pulse_shaping": True,
                "type_pulse_shape": "sine_squared",
                "time_pulse_shape_rolloff_us": 10.0,
                "freq_pulse_shape_sample_khz": 250.,


                "time_heating_us": 100.0,
                "att_phaser_db": 24.,
                "freq_global_offset_mhz": 0.0,
                "freq_osc_khz_list": "[-1284.697, 1284.697, 0.0, 0.0]",
                "ampl_osc_frac_list": "[5.94, 9.27, 11.89, 0.0]",
                "phase_osc_turns_list": "[0.3257, 0.0, 0., 0.0]",
                "phase_osc_ch1_offset_turns": "[0.0, 0.0, 0.5, 0.5, 0.5]",


                "enable_wav_scale": False,
                "target_wav_scale": "Time (Total)",
                "wav_osc_scale_list": {"sequence": [1.0], "ty": "ExplicitScan"},


                "enable_phase_shift_keying": False,
                "phase_osc0_psk_turns": "[0.0, 0.5]",
                "phase_osc1_psk_turns": "[0.0, 0.5]",
                "phase_osc2_psk_turns": "[0.0, 0.0]",
                "phase_osc3_psk_turns": "[0.0, 0.0]",
                "phase_osc4_psk_turns": "[0.0, 0.0]",


                "enable_psk_delay": False,
                "time_psk_delay_us_list": {
                    "sequence": [1000.0], "ty": "ExplicitScan"
                }
            },
        }

        return expid

    def get_scans(self):
        """
        fd
        """

        '''
        CONFIGURE
        '''
        base_arr = np.array([5, 7.8, 10, 0.])
        base_att_db = 31.5

        base_nbar_ref = 0.322
        base_time_us = 500.
        base_freq_span_khz = 8.
        base_freq_step_khz = 0.2

        fock_range = (0, 1, 2, 3)
        time_us_arr = [250]
        # time_us_arr = [4000, 5000, 6000, 7000, 8000, 9000]

        # target_nbar_arr = np.array([0.71136, 0.395536, 0.35568, 0.318154]) ** 2.
        target_nbar_arr = np.array([0.4, 0.4, 0.4, 0.4]) ** 2.


        '''
        PREPARE SCAN ARR
        '''
        # calculate ampl scaling array for each fock state
        ampl_scale_factor_arr = (target_nbar_arr / base_nbar_ref) ** (1. / 4.)

        # shuffle the time array
        shuffle(time_us_arr)


        '''
        CREATE SCAN ARR
        '''
        scan_config = []
        for time_us in time_us_arr:
            # get general scaling factor
            time_scale_factor = base_time_us / time_us

            # calculate new freq range & step
            freq_range_span = base_freq_span_khz * time_scale_factor
            freq_range_step = base_freq_step_khz * time_scale_factor

            # calculate new timings
            time_shape_rolloff_us = (base_time_us / 10.) / time_scale_factor

            # create temporary list so we can shuffle the fock_range
            dict_holder_list_tmp = []
            for fock_num in fock_range:
                # calculate new amplitude scale
                osc_ampl_arr = (base_arr * ampl_scale_factor_arr[fock_num]) * (time_scale_factor**0.5)

                # create expid update dict
                dict_dj = {
                    # general
                    'arguments.att_phaser_db': base_att_db,

                    # timing
                    'arguments.time_heating_us': time_us,
                    'arguments.time_pulse_shape_rolloff_us': time_shape_rolloff_us,

                    # freq
                    'arguments.freq_osc_sweep_khz_list.span': freq_range_span,
                    'arguments.freq_osc_sweep_khz_list.step': freq_range_step,

                    # fock-specific
                    'arguments.fock_num_prep': fock_num,
                    'arguments.fock_num_read': str(fock_num),
                    'arguments.ampl_osc_frac_list': str(list(osc_ampl_arr)),
                }

                dict_holder_list_tmp.append(dict_dj)

            # shuffle the fock ordering and update main scan_config
            shuffle(dict_holder_list_tmp)
            scan_config.extend(dict_holder_list_tmp)

        # # randomize the submitted configs before returning
        # shuffle(scan_config)
        return scan_config

    def run(self):
        # get expid dujour
        expid_target = self.get_expid()
        scan_list = self.get_scans()

        # batch submit scanned experiments
        for scan_target_params in scan_list:
            # create deep copy of expid
            expid_local = deepcopy(expid_target)

            # update expid with arguments to be scanned
            for param_key, param_val in scan_target_params.items():
                update_deep(expid_local, param_key, param_val)

            # submit expid to scheduler
            rid_dj = self.scheduler.submit(pipeline_name='main', expid=expid_local)


    '''
    STORAGE
    '''
    def _store(self):
        expid_laserscan = {
            "devarg_override": {},
            "log_level": 30,
            "repo_rev": "34150853bfa4366dc0a345ae93b40f0f7e3b60e8",
            "file": "experiments\\diagnostics\\LaserScan.py",
            "class_name": "LaserScan",

            "arguments": {
                "repetitions": 25,

                "enable_linetrigger": False,
                "time_linetrig_holdoff_ms_list": {
                    "sequence": [0.1], "ty": "ExplicitScan"
                },

                "freq_qubit_scan_mhz": {
                    "center": 100.3001, "span": 0.002, "step": 2e-05,
                    "randomize": True, "seed": None, "ty": "CenterScan"
                },
                "time_qubit_us": 3500.0,
                "ampl_qubit_pct": 30.0,
                "att_qubit_db": 31.5,
                "enable_pulseshaping": False,

                "enable_rescue": False,
                "enable_resuscitate": False,
                "enable_aperture": False,
                "add_397nm_spinpol": False,
                "repetitions_per_rescue": 1
            }
        }

        scan_list = (
            {
                'arguments.ampl_qubit_pct': 11,
                'arguments.freq_qubit_scan_mhz.center': 101.0856
            },


            {
                'arguments.ampl_qubit_pct': 19,
                'arguments.freq_qubit_scan_mhz.center': 100.7344
            },
            {
                'arguments.ampl_qubit_pct': 30,
                'arguments.freq_qubit_scan_mhz.center': 100.4437
            },
            {
                'arguments.ampl_qubit_pct': 30,
                'arguments.freq_qubit_scan_mhz.center': 100.2985
            },


            {
                'arguments.ampl_qubit_pct': 20,
                'arguments.freq_qubit_scan_mhz.span': 0.004,
                'arguments.freq_qubit_scan_mhz.center': 104.1370
            },
            {
                'arguments.ampl_qubit_pct': 30,
                'arguments.freq_qubit_scan_mhz.span': 0.004,
                'arguments.freq_qubit_scan_mhz.center': 110.2366
            },
            {
                'arguments.ampl_qubit_pct': 16,
                'arguments.freq_qubit_scan_mhz.span': 0.004,
                'arguments.freq_qubit_scan_mhz.center': 113.2865
            },
        )
