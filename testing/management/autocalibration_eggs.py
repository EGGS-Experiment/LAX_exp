import numpy as np
from artiq.experiment import *

from enum import Enum
from collections import deque
from asyncio import get_event_loop, Event

from sipyco import pyon
from sipyco.sync_struct import Subscriber
from sipyco.asyncio_tools import atexit_register_coroutine


class Status(Enum):
    experiment_submitting =     0
    experiment_waiting =        1
    calibration_initializing =  2
    calibration_waiting =       3
    calibration_processing =    4

def update_deep(dict_dj, keyspec, val):
    """
    Modify values within a nested dict using keys in dot form.
    Arguments:
        todo: document
    """
    keylist = keyspec.split('.')
    for key in keylist[:-1]:
        dict_dj = dict_dj[key]
    dict_dj[keylist[-1]] = val


class Autocalibration(EnvExperiment):
    """
    idk document i guess
    """

    """
    SETUP
    """
    def build(self):
        # get core devices
        self.ccb =              self.get_device("ccb")
        self.scheduler =        self.get_device("scheduler")

        # todo: maybe can make expid and stuff be uploaded from a file, and set the filename as an argument
        # autocalibration parameters
        self.experiments_per_calibration =  5
        self.experiment_repetitions =       150

    def prepare(self):
        # create necessary data structures
        self._status =                  Status.experiment_submitting
        self._running_experiments =     set()

        self._pending_calibrations =    deque()
        self._running_calibrations =    set()
        self._calibration_callback =    lambda x: x
        self._calibration_results =     dict()

        # prepare parameters and related values
        self._prepare_parameters()

        # create list of calibrations and experiments
        self._prepare_expids()

    def _prepare_parameters(self):
        """
        Prepare holder variables for autocalibration.
        """
        # create list of parameters for calibration (to prevent overriding of experiment parameters)
        self.calibration_parameters =       {
            # spectra
            'freq_carrier_mhz':             102.1931,
            'freq_rsb_laserscan_mhz':       101.4884,
            'freq_bsb_laserscan_mhz':       102.8958,

            # SBC
            'ampl_quench_pow_uw':           8.33,

            # readout
            'time_carrier_rabi_us':         11.4,
            'time_bsb_rabi_us':             128,

            # EGGS
            'freq_rsb_eggs_readout_mhz':    101.4911,
            'freq_bsb_eggs_readout_mhz':    102.8962,
            'freq_wsec_eggs_readout_khz':   1406.43
        }

        # create list of parameters to continually update the experiments with
        self.experiment_parameters =       {
            # sideband cooling
            'freq_sideband_cooling_mhz_pct_list':           pyon.encode({101.4873: 100.}),
            'ampl_quench_pct':                              2.25,

            # readout
            'freq_rsb_readout_mhz_list.sequence':           [101.4911],
            'freq_bsb_readout_mhz_list.sequence':           [102.8962],
            'time_readout_us_list.sequence':                [128],

            # EGGS
            'freq_eggs_heating_secular_khz_list.sequence':  [1406.43]
        }

    def _prepare_expids(self):
        """
        Prepare expids for calibration sequences and the main experiment queue.
        """
        '''
        Calibration List
        A deque that contains number of "calibration sets".
        Each calibration set runs the same expid for a set of parameters.
        
        parameter_name: the parameter of the calibration set that is to be adjusted.
        sweep_function: returns the parameters needed for the calibration.
                        e.g. returns the previous known carrier frequency if we want to run
                            a laser scan that finds the carrier.
        callback:       a function that is called when the entire calibration set is completed.
        '''
        # calibration set: dict(name, preprocess_func, postprocess_func, expid)
        self.calibrations_list =        deque([
            # calibration #1 - laser scan carrier
            {
                'name':             "ls_carrier",
                'preprocess_func':  self._ls_carrier_preprocess,
                'postprocess_func': self._ls_carrier_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\diagnostics\\LaserScan.py",
                    "class_name":   "LaserScan",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":          30,
                        "att_qubit_db":         31.5,
                        "ampl_qubit_pct":       38,
                        "time_qubit_us":        5000,
                        "freq_qubit_scan_mhz":  {"center": 102.1931, "span": 0.01, "step": 0.00005,
                                                 "randomize": True, "seed":None, "ty": "CenterScan"}
                    }
                }
            },
            # calibration #2 - laser scan RSB
            {
                'name':             "ls_rsb",
                'preprocess_func':  self._ls_rsb_preprocess,
                'postprocess_func': self._ls_rsb_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\diagnostics\\LaserScan.py",
                    "class_name":   "LaserScan",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":          20,
                        "att_qubit_db":         31.5,
                        "ampl_qubit_pct":       42,
                        "time_qubit_us":        5000,
                        "freq_qubit_scan_mhz":  {"center": 101.4881, "span": 0.01, "step": 0.00005,
                                                "randomize": True, "seed": None, "ty": "CenterScan"}
                    }
                }
            },
            # calibration #3 - laser scan BSB
            {
                'name':             "ls_bsb",
                'preprocess_func':  self._ls_bsb_preprocess,
                'postprocess_func': self._ls_bsb_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\diagnostics\\LaserScan.py",
                    "class_name":   "LaserScan",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":          20,
                        "att_qubit_db":         31.5,
                        "ampl_qubit_pct":       42,
                        "time_qubit_us":        5000,
                        "freq_qubit_scan_mhz":  {"center": 102.8975, "span": 0.01, "step": 0.00005,
                                                 "randomize": True, "seed": None, "ty": "CenterScan"}
                    }
                }
            },
            # calibration #4 - carrier rabi flop
            {
                'name':             "rabi_carrier",
                'preprocess_func':  self._rabi_carrier_preprocess,
                'postprocess_func': self._rabi_carrier_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\diagnostics\\RabiFlopping.py",
                    "class_name":   "RabiFlopping",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":          50,
                        "cooling_type":         "Doppler",
                        # readout
                        "time_rabi_us_list":    {"npoints": 200, "randomize": 2, "seed": None,
                                                 "start": 1.0, "stop": 50.0, "ty": "RangeScan"},
                        "freq_rabiflop_mhz":    102.1931,
                        "att_readout_db":       8.0,
                        # rescue
                        "rescue_enable":            False,
                        "repetitions_per_rescue":   1,
                        "resuscitate_ion":          False
                    }
                }
            },
            # calibration #5 - EGGS sideband readout
            {
                'name':             "eggs_sideband",
                'preprocess_func':  self._eggs_sideband_preprocess,
                'postprocess_func': self._eggs_sideband_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\eggs_heating\\EGGSHeating.py",
                    "class_name":   "EGGSHeating",
                    "log_level":    30,
                    "arguments": {
                        # config
                        "repetitions":      30,
                        "randomize_config": False,
                        "sub_repetitions":  1,
                        # SBC
                        "calibration_continuous":               False,
                        "sideband_cycles_continuous":           1,
                        "time_sideband_cooling_us":             10888.0,
                        "pct_per_spin_polarization":            40.0,
                        "freq_sideband_cooling_mhz_pct_list":   pyon.encode({101.4872: 100}),
                        "att_sidebandcooling_continuous_db":    8.0,
                        "ampl_quench_pct":                      2.08,
                        # readout
                        "freq_rsb_readout_mhz_list":    {"center": 101.4911, "span": 0.02, "step": 0.00025,
                                                         "randomize": True, "seed": None, "ty": "CenterScan"},
                        "freq_bsb_readout_mhz_list":    {"center": 102.8962, "span": 0.02, "step": 0.00025,
                                                         "randomize": True, "seed": None, "ty": "CenterScan"},
                        "ampl_sideband_readout_pct":    50.0,
                        "att_sideband_readout_db":      8.0,
                        "time_sideband_readout_us":     128.0,
                        # rescue
                        "rescue_enable":            False,
                        "repetitions_per_rescue":   1,
                        "resuscitate_ion":          False,
                        "add_397nm_spinpol":        False,
                        "death_detection":          True,
                        "time_readout_us_list":     {"sequence": [128.0], "ty": "ExplicitScan"},
                        # EGGS
                        "freq_eggs_heating_carrier_mhz_list":   {"sequence": [82.0], "ty": "ExplicitScan"},
                        "freq_eggs_heating_secular_khz_list":   {"sequence": [1406.66], "ty": "ExplicitScan"},
                        "enable_amplitude_calibration":         False,
                        "ampl_eggs_heating_rsb_pct":            60.0,
                        "ampl_eggs_heating_bsb_pct":            0.0,
                        "att_eggs_heating_db":                  3,
                        "time_eggs_heating_ms":                 1.0,
                        "phase_eggs_heating_rsb_turns_list":    {"sequence": [0.0], "ty": "ExplicitScan"},
                        "phase_eggs_heating_bsb_turns":         0.0,
                        "enable_pulse_shaping":                 False,
                        "type_pulse_shape":                     "sine_squared",
                        "time_pulse_shape_rolloff_us":          100.0,
                        "freq_pulse_shape_sample_khz":          500,
                        "enable_dynamical_decoupling":          True,
                        "ampl_eggs_dynamical_decoupling_pct":   0.05,
                        "enable_dd_phase_shift_keying":         False,
                        "num_dynamical_decoupling_phase_shifts": 3
                    }
                }
            },
            # calibration #6 - EGGS secular readout
            {
                'name':             "eggs_secular",
                'preprocess_func':  self._eggs_secular_preprocess,
                'postprocess_func': self._eggs_secular_postprocess,
                'expid': {
                    "file":         "LAX_exp\\experiments\\eggs_heating\\EGGSHeating.py",
                    "class_name":   "EGGSHeating",
                    "log_level":    30,
                    "arguments": {
                        # config
                        "repetitions":      30,
                        "randomize_config": False,
                        "sub_repetitions":  1,
                        # SBC
                        "calibration_continuous":               False,
                        "sideband_cycles_continuous":           1,
                        "time_sideband_cooling_us":             10888.0,
                        "pct_per_spin_polarization":            40.0,
                        "freq_sideband_cooling_mhz_pct_list":   pyon.encode({101.4872: 100}),
                        "att_sidebandcooling_continuous_db":    8.0,
                        "ampl_quench_pct":                      2.08,
                        # readout
                        "freq_rsb_readout_mhz_list":    {"sequence": [101.4911], "ty": "ExplicitScan"},
                        "freq_bsb_readout_mhz_list":    {"sequence": [102.8962], "ty": "ExplicitScan"},
                        "ampl_sideband_readout_pct":    50.0,
                        "att_sideband_readout_db":      8.0,
                        "time_sideband_readout_us":     128.0,
                        # rescue
                        "rescue_enable":            False,
                        "repetitions_per_rescue":   1,
                        "resuscitate_ion":          False,
                        "add_397nm_spinpol":        False,
                        "death_detection":          True,
                        "time_readout_us_list":     {"sequence": [128.0], "ty": "ExplicitScan"},
                        # EGGS
                        "freq_eggs_heating_carrier_mhz_list":   {"sequence": [82.0], "ty": "ExplicitScan"},
                        "freq_eggs_heating_secular_khz_list":   {"center": 1406.66, "span": 5.0, "step": 0.1,
                                                                 "randomize": True, "seed": None, "ty": "CenterScan"},
                        "enable_amplitude_calibration":         False,
                        "ampl_eggs_heating_rsb_pct":            60.0,
                        "ampl_eggs_heating_bsb_pct":            0.0,
                        "att_eggs_heating_db":                  3,
                        "time_eggs_heating_ms":                 1.0,
                        "phase_eggs_heating_rsb_turns_list":    {"sequence": [0.0], "ty": "ExplicitScan"},
                        "phase_eggs_heating_bsb_turns":         0.0,
                        "enable_pulse_shaping":                 False,
                        "type_pulse_shape":                     "sine_squared",
                        "time_pulse_shape_rolloff_us":          100.0,
                        "freq_pulse_shape_sample_khz":          500,
                        "enable_dynamical_decoupling":          True,
                        "ampl_eggs_dynamical_decoupling_pct":   0.05,
                        "enable_dd_phase_shift_keying":         False,
                        "num_dynamical_decoupling_phase_shifts": 3
                    }
                }
            }
        ])

        # create list of experiments to submit
        # note: each element in the deque should be a simple expid dict
        self.pending_experiments = deque([
        {
            "file":         "LAX_exp\\experiments\\eggs_heating\\EGGSHeating.py",
            "class_name":   "EGGSHeating",
            "log_level":    30,
            "arguments": {
                # config
                "repetitions":      30,
                "randomize_config": False,
                "sub_repetitions":  1,
                # SBC
                "calibration_continuous":               False,
                "sideband_cycles_continuous":           1,
                "time_sideband_cooling_us":             10888.0,
                "pct_per_spin_polarization":            40.0,
                "freq_sideband_cooling_mhz_pct_list":   pyon.encode({101.4872: 100}),
                "att_sidebandcooling_continuous_db":    8.0,
                "ampl_quench_pct":                      2.08,
                # readout
                "freq_rsb_readout_mhz_list":    {"sequence": [101.4911], "ty": "ExplicitScan"},
                "freq_bsb_readout_mhz_list":    {"sequence": [102.8962], "ty": "ExplicitScan"},
                "ampl_sideband_readout_pct":    50.0,
                "att_sideband_readout_db":      8.0,
                "time_sideband_readout_us":     128.0,
                # rescue
                "rescue_enable":            False,
                "repetitions_per_rescue":   1,
                "resuscitate_ion":          False,
                "add_397nm_spinpol":        False,
                "death_detection":          True,
                "time_readout_us_list":     {"sequence": [128.0], "ty": "ExplicitScan"},
                # EGGS
                "freq_eggs_heating_carrier_mhz_list":   {"center": 83.09, "span": 0.02, "step": 0.0002,
                                                         "randomize": True, "seed": None, "ty": "CenterScan"},
                "freq_eggs_heating_secular_khz_list":   {"sequence": [1406.66], "ty": "ExplicitScan"},
                "enable_amplitude_calibration":         False,
                "ampl_eggs_heating_rsb_pct":            60.0,
                "ampl_eggs_heating_bsb_pct":            0.0,
                "att_eggs_heating_db":                  3,
                "time_eggs_heating_ms":                 1.0,
                "phase_eggs_heating_rsb_turns_list":    {"sequence": [0.0], "ty": "ExplicitScan"},
                "phase_eggs_heating_bsb_turns":         0.0,
                "enable_pulse_shaping":                 False,
                "type_pulse_shape":                     "sine_squared",
                "time_pulse_shape_rolloff_us":          100.0,
                "freq_pulse_shape_sample_khz":          500,
                "enable_dynamical_decoupling":          True,
                "ampl_eggs_dynamical_decoupling_pct":   0.05,
                "enable_dd_phase_shift_keying":         False,
                "num_dynamical_decoupling_phase_shifts": 3
            }
            # "repo_rev": "41897f76abd02dd6cd37e01db49b9dc9fa4f111d"
        } for i in range(self.experiment_repetitions)])

        # shuffle experiment order
        np.random.shuffle(self.pending_experiments)


    """
    CALIBRATION PRE & POST-PROCESSING
    """
    def _ls_carrier_preprocess(self):
        return {"freq_qubit_scan_mhz.center": self.calibration_parameters['freq_carrier_mhz']}

    def _ls_carrier_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)

        # update calibration parameters with carrier value
        # note: result is list of [[peak_freq_mhz, peak_val_pop]]
        freq_carrier_mhz = results_dict['ls_carrier'][0, 0]
        self.calibration_parameters['freq_carrier_mhz'] = freq_carrier_mhz
        print('\t\tCarrier Calib')
        print('\t\t\tCarrier freq: {:.4f}'.format(freq_carrier_mhz))

    def _ls_rsb_preprocess(self):
        return {"freq_qubit_scan_mhz.center": self.calibration_parameters['freq_rsb_laserscan_mhz']}

    def _ls_rsb_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)
        print(results_dict)

        # update calibration_parameters with rsb value
        # note: result is list of [[peak_freq_mhz, peak_val_pop]]
        freq_rsb_mhz = results_dict['ls_rsb'][0, 0]
        self.calibration_parameters['freq_rsb_laserscan_mhz'] = freq_rsb_mhz

        # update SBC freq, RSB freq, and secular freq in experiment_parameters
        freq_sbc_mhz = freq_rsb_mhz - 0.0009
        self.experiment_parameters['freq_sideband_cooling_mhz_pct_list'] = pyon.encode({freq_sbc_mhz: 100.})

        freq_rsb_readout_mhz = freq_rsb_mhz + 0.0023
        self.experiment_parameters['freq_rsb_readout_mhz_list.sequence'] = [freq_rsb_readout_mhz]

        freq_wsec_khz = 2. * (self.calibration_parameters['freq_carrier_mhz'] - freq_rsb_mhz) * 1000. - 2.72
        self.experiment_parameters['freq_eggs_heating_secular_khz_list.sequence'] = [freq_wsec_khz]
        print('\t\tLS RSB Calib')
        print('\t\t\tSBC freq: {:.4f}\n\t\t\tRSB Readout: {:.4f}\n\t\t\twsec freq: {:.4f}'.format(freq_sbc_mhz,
                                                                                                  freq_rsb_readout_mhz,
                                                                                                  freq_wsec_khz))

    def _ls_bsb_preprocess(self):
        return {"freq_qubit_scan_mhz.center": self.calibration_parameters['freq_bsb_laserscan_mhz']}

    def _ls_bsb_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)

        # update calibration_parameters with bsb value
        # note: result is list of [[peak_freq_mhz, peak_val_pop]]
        freq_bsb_mhz = results_dict['ls_bsb'][0, 0]
        self.calibration_parameters['freq_bsb_laserscan_mhz'] = freq_bsb_mhz

        # update BSB readout freq in experiment_parameters
        freq_bsb_readout_mhz = freq_bsb_mhz - 0.0013
        self.experiment_parameters['freq_bsb_readout_mhz_list.sequence'] = [freq_bsb_readout_mhz]

        print('\t\tBSB Calib')
        print('\t\t\tBSB freq: {:.4f}'.format(freq_bsb_mhz))

    def _rabi_carrier_preprocess(self):
        return {"freq_rabiflop_mhz": self.calibration_parameters['freq_carrier_mhz']}

    def _rabi_carrier_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)
        # update calibration_parameters with bsb value
        # note: result is list of [[peak_freq_mhz, peak_val_pop]]
        time_readout_carrier_us = results_dict['rabi_carrier'][0][0]
        self.calibration_parameters['time_carrier_rabi_us'] = time_readout_carrier_us

        # update BSB readout freq in experiment_parameters
        time_readout_bsb_us = time_readout_carrier_us * 11.36
        self.experiment_parameters['time_readout_us_list.sequence'] = [time_readout_bsb_us]

        print('\t\tRabi Carrier Calib')
        print('\t\t\tCarrier rabi time: {:.4f}'.format(time_readout_carrier_us))

    def _eggs_sideband_preprocess(self):
        # get EGGS sideband values from experiment_parameters since they require processing
        freq_sbc_list =         self.experiment_parameters['freq_sideband_cooling_mhz_pct_list']
        time_readout_us =       self.experiment_parameters['time_readout_us_list.sequence']
        freq_eggs_rsb_mhz =     self.experiment_parameters['freq_rsb_readout_mhz_list.sequence']
        freq_eggs_bsb_mhz =     self.experiment_parameters['freq_bsb_readout_mhz_list.sequence']
        freq_eggs_wsec_khz =    self.experiment_parameters['freq_eggs_heating_secular_khz_list.sequence']
        return {"freq_sideband_cooling_mhz_pct_list":           freq_sbc_list,
                "time_readout_us.sequence":                     time_readout_us,
                "freq_rsb_readout_mhz_list.center":             freq_eggs_rsb_mhz,
                "freq_bsb_readout_mhz_list.center":             freq_eggs_bsb_mhz,
                "freq_eggs_heating_secular_khz_list.sequence":  freq_eggs_wsec_khz}

    def _eggs_sideband_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)
        # update calibration_parameters with parameters
        # note: result is list of [[n, n_err], [fit_params_rsb, fit_err_rsb], [fit_params_bsb, fit_err_bsb]]
        freq_rsb_eggs_mhz = results_dict['eggs_sideband'][1][0][1] / 2.
        freq_bsb_eggs_mhz = results_dict['eggs_sideband'][2][0][1] / 2.
        self.experiment_parameters['freq_rsb_readout_mhz_list.sequence'] = [freq_rsb_eggs_mhz]
        self.experiment_parameters['freq_bsb_readout_mhz_list.sequence'] = [freq_bsb_eggs_mhz]

        print('\t\tEGGS Sideband Calib')
        print('\t\t\tSideband Freqs: {:.4f}, {:.4f}'.format(freq_rsb_eggs_mhz, freq_bsb_eggs_mhz))

    def _eggs_secular_preprocess(self):
        # get EGGS secular values from experiment_parameters since they require processing
        freq_sbc_list =         self.experiment_parameters['freq_sideband_cooling_mhz_pct_list']
        time_readout_us =       self.experiment_parameters['time_readout_us_list.sequence']
        freq_eggs_rsb_mhz =     self.experiment_parameters['freq_rsb_readout_mhz_list.sequence']
        freq_eggs_bsb_mhz =     self.experiment_parameters['freq_bsb_readout_mhz_list.sequence']
        freq_eggs_wsec_khz =    self.experiment_parameters['freq_eggs_heating_secular_khz_list.sequence']
        return {"freq_sideband_cooling_mhz_pct_list":           freq_sbc_list,
                "time_readout_us.sequence":                     time_readout_us,
                "freq_rsb_readout_mhz_list.sequence":           freq_eggs_rsb_mhz,
                "freq_bsb_readout_mhz_list.sequence":           freq_eggs_bsb_mhz,
                "freq_eggs_heating_secular_khz_list.center":    freq_eggs_wsec_khz}

    def _eggs_secular_postprocess(self, results_dict):
        # todo: check for any errors (new value should be reasonably close to old)
        # update calibration_parameters with parameters
        # note: result is list of [[n, n_err], [fit_params_secular, fit_err_secular]]
        freq_eggs_secular_khz = results_dict['eggs_secular'][1][0][1] / 2.
        self.experiment_parameters['freq_eggs_heating_secular_khz_list.sequence'] = [freq_eggs_secular_khz]

        print('\t\tEGGS Secular Calib')
        print('\t\t\tSecular Freq.: {:.4f}'.format(freq_eggs_secular_khz))




    """
    MAIN SEQUENCE
    """
    def run(self):
        try:
            # set up target_builder function for scheduler client
            _scheduler_struct = dict()
            def _update_scheduler(x):
                _scheduler_struct.clear()
                _scheduler_struct.update(x)
                return _scheduler_struct

            # set up target_builder function for dataset client
            _dataset_struct = dict()
            def _update_datasets(x):
                _dataset_struct.clear()
                _dataset_struct.update(x)
                return _dataset_struct


            # set up event loop for subscribers
            loop =              get_event_loop()
            self.stop_event =   Event()

            # create subscribers
            self.scheduler_subscriber = Subscriber('schedule',
                                                   _update_scheduler,
                                                   lambda mod: self._process_scheduler_update(_scheduler_struct, mod))
            self.dataset_subscriber = Subscriber('datasets',
                                                 _update_datasets,
                                                 lambda mod: self._process_datasets_update(_dataset_struct, mod))
            # connect subscribers
            loop.run_until_complete(self.scheduler_subscriber.connect('::1', 3250))
            loop.run_until_complete(self.dataset_subscriber.connect('::1', 3250))

            # ensure subscribers close connection upon exit
            atexit_register_coroutine(self.scheduler_subscriber.close)
            atexit_register_coroutine(self.dataset_subscriber.close)

            # submit initial experiments first
            self._status = Status.experiment_submitting
            loop.call_later(2, self._submit_experiments)

            # run event loop indefinitely
            loop.run_until_complete(self.stop_event.wait())
            # loop.stop()
            # loop.close()

        except Exception as e:
            print('\t\t\tError during Run: {}'.format(repr(e)))
        finally:
            # loop.stop()
            # loop.close()
            print('\t-----------------------------AUTOCALIBRATION DONE-----------------------------')


    """
    SUBSCRIBER METHODS
    """
    def _process_scheduler_update(self, scheduler_dict, mod=None):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.

        Arguments:
            scheduler_dict  dict()  : ***todo: document***
        """
        # todo: make it work on ModActions instead of the actual dict
        # todo: try to minimize overhead by returning ASAP
        # see if any of our submitted experiments are running
        for exp_rid, exp_notification in scheduler_dict.items():
            # remove the experiment from our checklist if it has finished running
            if (exp_rid in self._running_experiments) and (exp_notification['status'] == 'deleting'):
                print('\t\tExperiment finished - RID: {:d}'.format(exp_rid))
                self._running_experiments.remove(exp_rid)

        # begin recalibration process if all submitted experiments have finished
        # and we have not already begun recalibrating
        if (len(self._running_experiments) == 0) and (self._status == Status.experiment_waiting):
            print('\tAutocalibration: EXPERIMENTS FINISHED ===> CALIBRATING')
            self._status = Status.calibration_initializing
            self._initialize_calibrations()

    def _process_datasets_update(self, dataset_dict, mod=None):
        """
        # todo: redocument
        Checks if any experiments are running and sends a Signal to clients accordingly.
        """
        # todo: check for ion death
        # todo: if ion death DURING a calibration, then rerun previous calibration
        # todo: can rerun by clearing values and restarting the stage
        # todo: if ion death DURING an experiment, then rerun previous experiment

        # ignore any updates not about completed calibrations
        try:
            mod_key = mod.get('key', [])
            # print(mod)
            if ('rid' not in mod_key) or (self._status != Status.calibration_waiting):
                return
        # needed in case mod['key'] is an int (which it is when the dataset/update is a list,
        # e.g. for mutate_dataset with temp.trace.counts
        except Exception as e:
            return

        # check if modification concerns one of our calibration experiments
        rid_num = mod['value'][1] if mod['action'] == 'setitem' else None
        if rid_num in self._running_calibrations:
            print('\t\tCalibration finished - RID: {:d}'.format(rid_num))

            # extract results from dataset_dict
            results_path = '.'.join(mod_key.split('.')[:-1] + ['results'])
            try:
                # store the results for processing by calibration later on
                self._calibration_results[rid_num]['results'] = dataset_dict[results_path][1]
            except KeyError as e:
                print("\t\tError during autocalib: unable to get results for RID: {:d}".format(rid_num))
            finally:
                # ensure we always remove rid key from _running_calibrations
                self._running_calibrations.remove(rid_num)

        # continue to calibration processing stage if all calibrations have finished
        if len(self._running_calibrations) == 0:
            print("\tAutocalibration: CALIBRATIONS FINISHED ===> PROCESSING")
            self._status = Status.calibration_processing
            self._process_calibration_stage()


    """
    EXPERIMENT STAGE
    """
    def _submit_experiments(self):
        """
        Submit a batch of experiments
        """
        # batch submit <num_exps> number of experiments at a time
        for i in range(self.experiments_per_calibration):

            # get expid from queue
            try:
                expid_dj = self.pending_experiments.popleft()
            except IndexError:
            # stop autocalibration if experiments are all completed
                self.stop_event.set()
                return

            # update expid with current parameters
            for param_key, param_val in self.experiment_parameters.items():
                update_deep(expid_dj['arguments'], param_key, param_val)

            # submit experiment to scheduler
            rid_dj = self.scheduler.submit(pipeline_name='calibrations', expid=expid_dj)
            self._running_experiments.update([rid_dj])
            print('\t\tSubmitted experiment - RID: {:d}'.format(rid_dj))

        # change status to waiting
        self._status = Status.experiment_waiting


    """
    CALIBRATION STAGE
    """
    def _initialize_calibrations(self):
        """
        todo: document
        """
        print('\tAutocalibration: CALIBRATIONS INITIALIZING')
        # create a copy of the calibrations to be submitted
        self._pending_calibrations = self.calibrations_list.copy()
        self._submit_calibration_stage()

    def _submit_calibration_stage(self):
        """
        Submit a set of calibrations (a calibration stage).
        todo: document
        """
        '''
        PREPARE
        '''
        print('\tAutocalibration: CALIBRATIONS SUBMITTING')
        
        # todo: add error handling to check that pending_calibrations is non-empty
        # clear loop iterators and get new calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        calibration_stage = self._pending_calibrations.popleft()

        # extract values to set up this calibration stage
        calibration_name =              calibration_stage['name']
        _calibration_preprocess =       calibration_stage['preprocess_func']
        self._calibration_postprocess = calibration_stage['postprocess_func']

        try:

            # get calibraton expid and update it deeply with target parameters
            # (_calibration_preprocess should return dict where keys are expid parameter names and values are values)
            expid_dj = calibration_stage['expid'].copy()
            for parameter_name, parameter_value in _calibration_preprocess().items():
                update_deep(expid_dj['arguments'], parameter_name, parameter_value)
        except Exception as e:
            print("\t\t\tproblem: {}".format(repr(e)))
            print('\t\t\t\tparam name: {}'.format(parameter_name))
            print('\t\t\t\tparam val: {}'.format(parameter_value))


        '''
        SUBMIT CALIBRATIONS
        '''
        # submit calibrated expid to scheduler and update holding structures
        rid_dj = self.scheduler.submit(pipeline_name='calibrations', expid=expid_dj)
        self._running_calibrations.update([rid_dj])
        self._calibration_results[rid_dj] = {
            'name':     calibration_name,
            'results':  None
        }
        print('\t\tSubmitting calibration - RID: {:d}'.format(rid_dj))

        # change status to calibration_waiting
        self._status = Status.calibration_waiting

    def _process_calibration_stage(self):
        """
        todo: document
        """
        # process _calibration_results to use the calibration name as the key (instead of its rid)
        # todo: add error handling if there's some problem with the results (e.g. bad fit to results)
        _calibration_results_sorted =   {result_dict['name']: result_dict['results']
                                         for exp_num, result_dict in self._calibration_results.items()}
        # pass calibration results to the associated callback for processing and updating
        self._calibration_postprocess(_calibration_results_sorted)


        # clean up calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        self._calibration_postprocess = None


        # check if we have more calibration sets to do, then load them into active calibration queue and submit again
        if len(self._pending_calibrations) != 0:
            print('\tAutocalibration: CALIBRATION PROCESSING FINISHED ===> NEXT CALIBRATION STAGE')
            self._submit_calibration_stage()
        # otherwise, change status and return to experiment stage
        else:
            self._status = Status.experiment_submitting
            print('\tAutocalibration: CALIBRATIONS FINISHED ===> SUBMITTING EXPERIMENTS')
            self._submit_experiments()


    """
    EXPID STORAGE
    """
    def __expid_storage(self):
        self._sbc_expids = deque([{
            "log_level": 30,
            "file": "LAX_exp\\experiments\\SidebandCooling.py",
            "class_name": "SidebandCooling",
            "arguments": {
                "repetitions": 20,
                "cooling_type": "Continuous",
                "freq_rsb_scan_mhz": {"center": 102.645, "randomize": True, "seed": None, "span": 0.02,
                                        "step": 0.0005, "ty": "CenterScan"},
                "freq_bsb_scan_mhz": {"center": 103.730, "randomize": True, "seed": None, "span": 0.02,
                                        "step": 0.0005, "ty": "CenterScan"},
                "time_readout_pipulse_us": 120.0,
                "ampl_readout_pipulse_pct": 50.0,
                "att_readout_db": 8.0,
                "calibration_continuous": False,
                "sideband_cycles_continuous": 1,
                "time_sideband_cooling_us": 8000.0,
                "pct_per_spin_polarization": 20.0,
                "freq_sideband_cooling_mhz_pct_list": "{102.647: 100}",
                "att_sidebandcooling_continuous_db": 8.0,
                "ampl_quench_pct": 4.0,
                "rescue_enable": False
            }
        } for i in range(self.experiment_repetitions)])
        
