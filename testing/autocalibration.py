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
    """
    keylist = keyspec.split('.')
    for key in keylist[:-1]:
        dict_dj = dict_dj[key]
    dict_dj[keylist[-1]] = val


class Autocalibration(EnvExperiment):
    """
    idk document i guess
    """

    def build(self):
        # get core devices
        self.ccb =              self.get_device("ccb")
        self.scheduler =        self.get_device("scheduler")

        # todo: maybe can make expid and stuff be uploaded from a file, and set the filename as an argument
        # autocalibration parameters
        self.experiments_per_calibration =  1
        self.experiment_repetitions =       2


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
        # tmp remove
        # set intermediate values for experiment
        self._freq_carrier_aom_mhz =    103.1885
        self._freq_axial_abs_mhz =      0.682
        self._freq_rf_abs_mhz =         0.788
        self._freq_eggs_abs_mhz =       1.086
        # tmp remove

        # create list of parameters for calibration (to prevent overriding of experiment parameters)
        self.calibration_parameters =       {
            'freq_qubit_scan_mhz.center':   self._freq_carrier_aom_mhz
        }

        # create list of parameters to continually update the experiments with
        self.experiment_parameters =       {
            'freq_rabiflop_mhz':                    self._freq_carrier_aom_mhz,
            'freq_sideband_cooling_mhz_pct_list':   pyon.encode({102.8485: 25, 102.793: 40, 102.6425: 35})
        }

    def _prepare_expids(self):
        # create list of calibrations to submit
        # need: expid, parameter_name, sweep_function, callback
        self.calibrations_list =        deque([
            {
                'parameter_name':   'freq_qubit_scan_mhz.center',
                'sweep_function':   self.sweep_func_1,
                'callback':         self.process_func_1,
                'expid': {
                    "file": "LAX_exp\\experiments\\LaserScan.py",
                    # "file":         "LAX_exp\\testing\\_autocalib_ls_test.py",
                    "class_name": "LaserScan",
                    # "class_name":   "autocalib_ls_test",
                    "log_level":    30,
                    "arguments": {
                        "repetitions":  25,
                        "att_qubit_db": 30.5,
                        "freq_qubit_scan_mhz": {
                            "center":       103.1880,
                            "span":         0.015,
                            "step":         0.0005,
                            "randomize":    True,
                            "seed":         None,
                            "ty":           "CenterScan"
                        }
                    }
                }
            }
            # {
            #     'parameter_name':   'freq_rabiflop_mhz',
            #     'sweep_function':   self.sweep_func_2,
            #     'callback':         self.process_func_2,
            #     'expid':{
            #         "log_level": 30,
            #         # "file": "experiments\\RabiFlopping.py",
            #         "file":         "LAX_exp\\testing\\_autocalib_rabiflop_test.py",
            #         # "class_name": "RabiFlopping",
            #         "class_name": "autocalib_rabiflop_test",
            #         "arguments": {
            #             "repetitions": 50,
            #             "cooling_type": "Doppler",
            #             "time_rabi_us_list": {
            #                 "npoints": 100,
            #                 "randomize": 2,
            #                 "seed": None,
            #                 "start": 1.0,
            #                 "stop": 40.0,
            #                 "ty": "RangeScan"
            #             },
            #             "freq_rabiflop_mhz": 103.2255,
            #             "att_readout_db": 8.0,
            #             "calibration_continuous": False,
            #             "sideband_cycles_continuous": 20,
            #             "time_sideband_cooling_us": 36000.0,
            #             "pct_per_spin_polarization": 20.0,
            #             "freq_sideband_cooling_mhz_pct_list": "{102.85: 20, 102.422: 40, 102.338: 40}",
            #             "att_sidebandcooling_continuous_db": 8.0,
            #             "ampl_quench_pct": 5.0,
            #             "rescue_enable": False, "repetitions_per_rescue": 1
            #         }
            #     }
            # }
        ])

        # create list of experiments to submit
        # note: each element in the deque should be a simple expid dict
        self.pending_experiments = deque([{
            "log_level":    30,
            "file":         "LAX_exp\\experiments\\RabiFlopping.py",
            "class_name":   "RabiFlopping",
            "arguments": {
                "repetitions":  25,
                "cooling_type": "SBC - Continuous",
                "time_rabi_us_list": {"npoints": 200, "randomize": 2, "seed": None, "start": 1.0, "stop": 200.0,
                                       "ty": "RangeScan"},
                "freq_rabiflop_mhz": 103.1885,
                "att_readout_db": 8.0,
                "calibration_continuous": False,
                "sideband_cycles_continuous": 20,
                "time_sideband_cooling_us": 35000.0,
                "pct_per_spin_polarization": 20.0,
                "freq_sideband_cooling_mhz_pct_list": pyon.encode({102.8485: 25, 102.793: 40, 102.6425: 35}),
                "att_sidebandcooling_continuous_db": 8.0,
                "ampl_quench_pct": 4.0,
                "rescue_enable": False
            }
        } for i in range(self.experiment_repetitions)])
        np.random.shuffle(self.pending_experiments)


    '''
    Calibration Processing
    '''

    def sweep_func_1(self, parameter_current):
        return [parameter_current,
                parameter_current - self._freq_axial_abs_mhz/2.,
                parameter_current - self._freq_rf_abs_mhz/2.,
                parameter_current - self._freq_eggs_abs_mhz/2.]

    def process_func_1(self, results_list):
        # guess new frequencies as mean of returned peak frequencies
        new_peak_values = np.array([np.mean(results[1][:, 0]) for results in results_list])
        print('new peak values: {}'.format(new_peak_values))

        # update values for calibration-only purposes
        self._freq_carrier_aom_mhz = new_peak_values[0]
        self._freq_axial_abs_mhz, self._freq_rf_abs_mhz, self._freq_eggs_abs_mhz = 2. * (new_peak_values[0] - new_peak_values[1:])

        # update calibration and experimental parameters
        self.calibration_parameters['freq_qubit_scan_mhz.center'] =         new_peak_values[0]
        self.experiment_parameters['freq_rabiflop_mhz'] =                   new_peak_values[0]
        self.experiment_parameters['freq_sideband_cooling_mhz_pct_list'] =  pyon.encode({new_peak_values[1] - 0.0007: 25.,
                                                                                         new_peak_values[2] - 0.00025: 40,
                                                                                         new_peak_values[3] - 0.00125: 35})

    def sweep_func_2(self, parameter_current):
        return np.linspace(parameter_current-0.02, parameter_current+0.02, 3)

    def process_func_2(self, results_list):
        print('\t\tresults 2: {}'.format(results_list))
        # todo: update self.experiment_parameters with freq_sideband_cooling_mhz_pct_list in pyon form


    '''
    MAIN SEQUENCE
    '''

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
            print('\t\t\tError during Run: {}'.format(e))
        finally:
            # loop.stop()
            # loop.close()
            print('\t-----------------------------AUTOCALIBRATION DONE-----------------------------')


    '''
    Subscriber Methods
    '''

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
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        # ignore any updates not about completed calibrations
        mod_key = mod.get('key', [])
        if ('rid' not in mod_key) or (self._status != Status.calibration_waiting):
            return

        # check if modification concerns one of our calibration experiments
        rid_num = mod['value'][1] if mod['action'] == 'setitem' else None
        if rid_num in self._running_calibrations:
            print('\t\tCalibration finished - RID: {:d}'.format(rid_num))

            # extract results from dataset_dict
            results_path = '.'.join(mod_key.split('.')[:-1] + ['results'])
            try:
                self._calibration_results[rid_num]['results'] = dataset_dict[results_path][1]
            except KeyError as e:
                print("\t\tError during autocalib: unable to get results for RID: {:d}".format(rid_num))
            finally:
                # ensure we always remove rid key from _running_calibrations
                self._running_calibrations.remove(rid_num)

        # continue to calibration processing stage if all calibrations have finished
        if len(self._running_calibrations) == 0:
            # print("\tAutocalibration: CALIBRATIONS FINISHED ===> PROCESSING")
            self._status = Status.calibration_processing
            self._process_calibration_stage()


    '''
    Experiment Methods
    '''

    def _submit_experiments(self):
        """
        todo: document
        """
        # batch submit <num_exps> number of experiments at a time
        for i in range(self.experiments_per_calibration):

            # get expid from queue
            try:
                expid_dj = self.pending_experiments.popleft()
            except IndexError:
                self.stop_event.set()
                return

            # update expid with current parameters
            [update_deep(expid_dj['arguments'], key_param, val_param)
             for key_param, val_param in self.experiment_parameters.items()]

            # submit experiment to scheduler
            rid_dj = self.scheduler.submit(pipeline_name='calibrations', expid=expid_dj)
            self._running_experiments.update([rid_dj])

            print('\t\tSubmitted experiment - RID: {:d}'.format(rid_dj))

        # change status to waiting
        self._status = Status.experiment_waiting


    '''
    Calibration Methods
    '''

    def _initialize_calibrations(self):
        """
        todo: document
        """
        # print('\tAutocalibration: CALIBRATIONS INITIALIZING')
        self._pending_calibrations = self.calibrations_list.copy()
        self._submit_calibration_stage()

    def _submit_calibration_stage(self):
        """
        todo: document
        """
        # print('\tAutocalibration: CALIBRATIONS SUBMITTING')
        
        # todo: add error handling to check that pending_calibrations is non-empty
        # clear loop iterators and get new calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        calibration_stage = self._pending_calibrations.popleft()

        # extract values to set up this calibration stage
        parameter_name =                calibration_stage['parameter_name']
        self._calibration_callback =    calibration_stage['callback']
        param_sweep_func =              calibration_stage['sweep_function']
        parameter_current_value =       self.calibration_parameters[parameter_name]

        # submit experiments to scan parameter around previously calibrated value
        for parameter_test_value in param_sweep_func(parameter_current_value):

            # get expid and update it deeply with the calibration test value
            expid_dj = calibration_stage['expid'].copy()
            update_deep(expid_dj['arguments'], parameter_name, parameter_test_value)

            # submit calibrated expid to scheduler and update holding structures
            rid_dj = self.scheduler.submit(pipeline_name='calibrations', expid=expid_dj)
            self._running_calibrations.update([rid_dj])
            self._calibration_results[rid_dj] = {
                'parameter_value':  parameter_test_value,
                'results':          None
            }
            print('\t\tSubmitting calibration - RID: {:d}'.format(rid_dj))

        # change status to waiting
        self._status = Status.calibration_waiting

    def _process_calibration_stage(self):
        """
        todo: document
        """
        # create a 2D array of [param_val, res] to pass to callback
        # todo: add error handling if there's some problem with the results
        _calibration_results = [(result_dict['parameter_value'], result_dict['results'])
                                for result_dict in self._calibration_results.values()]
        self._calibration_callback(_calibration_results)

        # clean up calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        self._calibration_callback = lambda x: x

        # check if we have more calibration sets to do, then load them into active calibration queue and submit again
        # otherwise, change status
        if len(self._pending_calibrations) == 0:
            self._status = Status.experiment_submitting
            print('\tAutocalibration: CALIBRATIONS FINISHED ===> SUBMITTING EXPERIMENTS')
            self._submit_experiments()
        else:
            # print('\tAutocalibration: CALIBRATION PROCESSING FINISHED ===> NEXT CALIBRATION STAGE')
            self._submit_calibration_stage()

    '''
    EXPID Storage
    '''

    def __expid_storage(self):
        self._eggsheating_expids = deque([{
            "log_level": 30,
            "file": "LAX_exp\\experiments\\EGGSHeating.py",
            "class_name": "EGGSHeating",
            "arguments": {
                "repetitions": 50,
                "cooling_type": "Continuous",
                "freq_rsb_scan_mhz": {"center": 102.6572, "randomize": True, "seed": None,
                                      "span": 0.01, "step": 0.0005, "ty": "CenterScan"},
                "freq_bsb_scan_mhz": {"center": 103.7389, "randomize": True, "seed": None,
                                      "span": 0.01, "step": 0.0005, "ty": "CenterScan"},
                "time_readout_pipulse_us": 110.0,
                "ampl_readout_pipulse_pct": 50.0,
                "att_readout_db": 8.0,
                "calibration_continuous": False,
                "sideband_cycles_continuous": 1,
                "time_sideband_cooling_us": 12000.0,
                "pct_per_spin_polarization": 20.0,
                "freq_sideband_cooling_mhz_pct_list": "{102.6572: 100}",
                "att_sidebandcooling_continuous_db": 8.0,
                "ampl_quench_pct": 4.0,
                "randomize_config": True,
                "freq_eggs_heating_carrier_mhz_list": {"sequence": [82.0], "ty": "ExplicitScan"},
                "freq_eggs_heating_secular_khz_list": {"center": 1087.7, "randomize": True, "seed": None,
                                                       "span": 5.0, "step": 0.25, "ty": "CenterScan"},
                "enable_amplitude_calibration": False,
                "ampl_eggs_heating_rsb_pct": 40.0,
                "ampl_eggs_heating_bsb_pct": 40.0,
                "att_eggs_heating_db": 5.0,
                "time_eggs_heating_ms": 1.0,
                "phase_eggs_heating_rsb_turns": 0.33,
                "phase_eggs_heating_bsb_turns": 0.25,
                "enable_pulse_shaping": False,
                "enable_dynamical_decoupling": True,
                "ampl_eggs_dynamical_decoupling_pct": 0.35,
                "enable_dd_phase_shift_keying": False,
                "num_dynamical_decoupling_phase_shifts": 3,
                "enable_dd_active_cancel": False
            }
        } for idk in range(5)])

        self._rabiflopping_expids = deque([{
            "log_level": 30,
            "file": "LAX_exp\\experiments\\RabiFlopping.py",
            "class_name": "RabiFlopping",
            "arguments": {
                "repetitions": 30,
                "cooling_type": "SBC - Continuous",
                "time_rabi_us_list": {"npoints": 100, "randomize": 2, "seed": None, "start": 1.0, "stop": 100.0,
                                       "ty": "RangeScan"},
                "freq_rabiflop_mhz": 103.1885,
                "att_readout_db": 8.0,
                "calibration_continuous": False,
                "sideband_cycles_continuous": 20,
                "time_sideband_cooling_us": 35000.0,
                "pct_per_spin_polarization": 20.0,
                "freq_sideband_cooling_mhz_pct_list": pyon.encode({102.8485: 25, 102.793: 40, 102.6425: 35}),
                "att_sidebandcooling_continuous_db": 8.0,
                "ampl_quench_pct": 4.0,
                "rescue_enable": False
            }
        } for idk in range(5)])

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
        
