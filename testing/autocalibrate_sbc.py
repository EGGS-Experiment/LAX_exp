import numpy as np
from artiq.experiment import *

from enum import Enum
from time import sleep
from queue import Queue
from collections import deque
from sipyco.sync_struct import Subscriber
from asyncio import get_event_loop, set_event_loop, Event
from sipyco.asyncio_tools import atexit_register_coroutine

# todo: support for multiple layers of calibrations (e.g. laser scan, then quench, then sideband cooling)
# todo: think about using pyonvalue instead of the strange rid abomination we have now


class Status(Enum):
    experiment_submitting =     0
    experiment_waiting =        1
    calibration_initializing =  2
    calibration_waiting =       3
    calibration_processing =    4


class Autocalibration(EnvExperiment):
    def build(self):
        # get core devices
        self.ccb =              self.get_device("ccb")
        self.scheduler =        self.get_device("scheduler")
        self.setattr_argument("submit_length",                           Scannable(
                                                                            default=CenterScan(40, 80, 1, randomize=True),
                                                                            global_min=0, global_max=300, global_step=1,
                                                                            unit="V", scale=1, ndecimals=3
                                                                        ))

        # autocalibration parameters
        self.num_exps =             3

        # base params - laser scan
        self.freq_carrier_mhz =     103.202
        self.freq_secular_khz =     1088.0
        self.freq_ls_scan_khz =     25.0

        # base params - sbc
        self.sbc_freq_span_khz =    4.
        self.sbc_freq_points =      5.



    def prepare(self):
        # create necessary data structures
        self.current_parameters =       dict()
        self._status =                  Status.experiment_waiting

        self._pending_experiments =     deque()
        self._running_experiments =     set()

        self._calibration_stage =       0
        self._pending_calibrations =    deque()
        self._running_calibrations =    set()
        self._calibration_results =     dict()

        # todo: idk
        self._calibration_callback =    lambda x: x
        # tmp remove

        # tmp remove
        # self._

        # create calibration expid list
        self.calibrations_list =        deque()
        # need: expid, parameter_name, sweep_function, callback

        # self.calibration_list = [
        #     {
        #         'status': 'pending',
        #         'callback':
        #         'expid': {
        #             "file": "LAX_exp\\experiments\\LaserScan.py",
        #             "class_name": "LaserScan",
        #             "log_level": 30,
        #             "arguments": {
        #                 "repetitions": 20,
        #                 "freq_qubit_scan_mhz": {
        #                     "center": 103.210,
        #                     "span": 0.02,
        #                     "step": 0.0005,
        #                     "randomize": True,
        #                     "seed": None,
        #                     "ty": "CenterScan"
        #                 }
        #             }
        #         }
        #     },
        #     {
        #
        #     }
        #
        # ]

        # create exp list
        self.pending_expid_list = [{
                "file":             "LAX_exp\\experiments\\LaserScan.py",
                "class_name":       "LaserScan",
                "log_level":        30,
                "arguments": {
                    "repetitions":      20,
                    "freq_qubit_scan_mhz": {
                        "center":       103.210,
                        "span":         0.02,
                        "step":         0.0005,
                        "randomize":    True,
                        "seed":         None,
                        "ty":           "CenterScan"
                    }
                }
            } for scan_val_v in self.scan_range_v]


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
            loop = get_event_loop()
            stop_event = Event()

            # create subscribers
            self.scheduler_subscriber = Subscriber('schedule',
                                                   _update_scheduler,
                                                   lambda mod: self._process_scheduler_update(_scheduler_struct, mod))
            # todo: ensure mod has nothing to do with datasets
            self.dataset_subscriber = Subscriber('datasets',
                                                 _update_datasets,
                                                 lambda mod: self._process_datasets_update(_dataset_struct, mod))
            # connect subscribers
            loop.run_until_complete(self.scheduler_subscriber.connect('::1', 3250))
            loop.run_until_complete(self.dataset_subscriber.connect('::1', 3250))

            # ensure subscribers close connection upon exit
            atexit_register_coroutine(self.scheduler_subscriber.close)
            atexit_register_coroutine(self.dataset_subscriber.close)

            # todo: submit initial experiments first

            # run event loop indefinitely
            loop.run_until_complete(stop_event.wait())
            loop.close()

        except Exception as e:
            print('\n\t\tError: {}'.format(e))
            raise


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
                self._running_experiments.remove(exp_rid)

        # begin recalibration process if all submitted experiments have finished
        # and we have not already begun recalibrating
        if (len(self._running_experiments) == 0) and (self._status == Status.experiment_waiting):
            print('\t\t\tAutocalibration:\tEXPERIMENTS FINISHED ===> CALIBRATING')
            self._status = Status.calibration_initializing
            self._initialize_calibrations()

    def _process_datasets_update(self, dataset_dict, mod=None):
        """
        # todo: redocument
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        mod_path = mod.get('path', [])

        # ignore any updates not about completed calibrations
        if ('rid' not in mod_path) or (self._status != Status.calibration_waiting):
            return

        rid_num = None
        if mod['actions'] == 'setitem':
            rid_num = mod['value'][1]

        # check if modification concerns one of our calibration experiments
        if rid_num in self._running_calibrations:

            # extract results from dataset_dict
            results_path = '.'.join(mod_path[:-1] + ['results'])
            try:
                calib_res = dataset_dict[results_path]
                self._calibration_results[rid_num]['results'] = calib_res
            except KeyError:
                print("\t\t\t\tError during autocalib: unable to get results for RID: {:d}".format(rid_num))
            finally:
                # ensure we always remove rid key from _running calibrations
                # to prevent forever loops
                self._running_calibrations.remove(rid_num)

        # continue to calibration processing stage if all calibrations have finished
        if len(self._running_calibrations) == 0:
            print("\t\t\tAutocalibration:\tCALIBRATIONS FINISHED ===> PROCESSING")
            self._status = Status.calibration_processing
            self._process_calibrations()


    '''
    Experiment Methods
    '''

    def _submit_experiments(self):
        """
        todo: document
        """
        # batch submit <num_exps> number of experiments at a time
        for i in range(self.num_exps):

            # get expid from queue with current calibrated parameters
            expid_dj = self.pending_experiments.popleft()
            expid_dj.update(self.current_parameters)

            # submit experiment to scheduler
            rid_dj = self._scheduler.submit(pipeline_name='main', expid=expid_dj)
            self._running_experiments.update(rid_dj)

        # change status to waiting
        self._status = Status.experiment_waiting


    '''
    Calibration Methods
    '''

    def _initialize_calibrations(self):
        """
        todo: document
        """
        self._pending_calibrations = self.calibrations_list.copy()
        self._submit_calibration_stage()


    def _submit_calibration_stage(self):
        """
        todo: document
        """
        # todo: add error handling to check that pending_calibrations is non-empty
        # clear loop iterators and get new calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        calibration_stage = self._pending_calibrations.popleft()

        # todo: document this section better, comments are off-ish
        # update callback method to be called by _process_calibration_stage
        # get current value of calibration parameter and
        # generate a parameter scan around previously calibrated value
        self._calibration_callback = calibration_stage['callback']
        parameter_name = self.calibrated_parameters[calibration_stage['parameter_name']]
        param_sweep_func = calibration_stage['scan_function']

        # create all expids for current calibration stage
        expid_list = {
            parameter_value: dict(calibration_stage['expid'],
                                  **{parameter_name: parameter_value})
            for parameter_value in param_sweep_func(parameter_name)
        }

        # submit calibrated expid to scheduler and update holding structures
        for parameter_value, expid in expid_list:
            rid_dj = self._scheduler.submit(pipeline_name='main', expid=expid_dj)
            self._running_calibrations.update(rid_dj)
            self._calibration_results.update({
                'rid':              rid_dj,
                'parameter_value':  parameter_value,
                'results':          None
            })

        # change status to waiting
        self._status = Status.calibration_waiting

    def _process_calibration_stage(self):
        """
        todo: document
        """
        # create a 2D array of [param_val, res] to pass to callback
        _calibration_results = [(result_dict['parameter_value'], result_dict['results'])
                                for result_dict in self._calibration_results]
        self._calibration_callback(_calibration_results)

        # clean up calibration stage
        self._running_calibrations.clear()
        self._calibration_results.clear()
        self._calibration_callback = lambda x: x

        # check if we have more calibration sets to do, then load them into active calibration queue and submit again
        # otherwise, change status
        if len(self._pending_calibrations) == 0:
            self._status = Status.experiment_submitting
            self._submit_experiments()
        else:
            self._submit_calibration_stage()

    '''
    User-subclassable Functions
    '''

    def generate_parameter_scan(self, parameter_dict):
        """
        todo: document
        """
        pass
