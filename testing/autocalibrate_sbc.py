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


class Status(Enum):
    waiting =       0
    calibrating =   1
    processing =    2


class Autocalibration(EnvExperiment):
    def build(self):
        self.ccb =              self.get_device("ccb")
        self.scheduler =        self.get_device("scheduler")
        self.setattr_argument("scan_range_v",                           Scannable(
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
        self._status =                  Status.waiting
        self.pending_expid_list =       deque()
        self.submitted_rid_list =       set()

        # tmp remove
        # data structure to hold calibrated parameters
        self._calibrated_params =       dict()
        self._calib_dataset_search =    dict()
        # tmp remove

        # create necessary values
        self.sbc_freq_list_mhz = np.ones(3) * (self.freq_carrier_mhz - self.freq_secular_khz / 2000.)

        # create exp list
        self.scan_range_v = list(self.scan_range_v)
        self.expid_list = [{
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
                                                   lambda mod: self._process_scheduler_update(_scheduler_struct))
            self.dataset_subscriber = Subscriber('datasets',
                                                 _update_datasets,
                                                 lambda mod: self._process_datasets_update(_dataset_struct))
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

    def _process_scheduler_update(self, scheduler_dict):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.

        Arguments:
            scheduler_dict  dict()  : ***todo: document***
        """
        # see if any of our submitted experiments are running
        for exp_rid, exp_notification in scheduler_dict.items():

            # remove the experiment from our checklist if it has finished running
            if (exp_rid in self.submitted_rid_list) and (exp_notification['status'] == 'deleting'):
                self.submitted_rid_list.remove(exp_rid)

        # begin recalibration process if all submitted experiments have finished
        # and we have not already begun recalibrating
        if (len(self.submitted_rid_list) == 0) and (self._status == Status.waiting):
            print('\t\t\tAutocalibration:\tEXPERIMENTS FINISHED ===> CALIBRATING')
            self._status = Status.calibrating
            self._submit_calibrations()

    def _process_datasets_update(self, dataset_dict):
        """
        # todo: redocument
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        # see if any of our pending calibrations have sent their results to the dataset manager
        for calib_rid, calib_dict in self._calib_dataset_search.items():

            # only check results of pending calibrations
            if calib_dict['status'] == 'Completed':
                continue

            # get keys to retrieve calibration results from dataset manager
            res_rid_key, res_result_key = ('{:s}.rid'.format(calib_dict['key']),
                                           '{:s}.results'.format(calib_dict['key']))

            # retrieve and store results if they are from one a calibration experiment we submitted
            if dataset_dict.get(res_rid_key, None) == calib_rid:
                calib_dict['results'] = dataset_dict.get(res_result_key)
                calib_dict['status'] = 'Completed'


        # get AND of completion status of all calibrations
        _calibration_completed_flag = all([
            calib_dict['status'] == 'Completed'
            for calib_dict in self._calib_dataset_search.values()
        ])

        # continue to calibration processing stage if all calibrations have finished
        if (self._status == Status.calibrating) and (_calibration_completed_flag is True):
            print('\t\t\tAutocalibration:\tCALIBRATIONS FINISHED ===> PROCESSING')
            self._status = Status.processing
            self._process_calibrations()


    '''
    Autocalibration Methods
    '''

    def _submit_experiments(self):
        """
        todo: document
        """
        # batch submit <num_exps> number of experiments at a time
        for i in range(self.num_exps):

            # get experiment expid from experiment queue
            # and update with our new calibrated parameters
            expid_dj = self.pending_expid_list.popleft()
            expid_dj.update(self._calibrated_params)

            # submit calibrated expid to scheduler and update
            rid_dj = self._scheduler.submit(pipeline_name='main', expid=expid_dj)
            self._calib_dataset_search[rid_dj] = {
                # todo: ensure we have the correct key for the experiment
                'key':          'idk',
                'status':       'pending',
                'results':      dict()
            }

        # change status to waiting
        self._status = Status.waiting

    def _submit_calibrations(self):
        """
        todo: document
        """
        # todo: get previous calibrated values
        # todo: scan a given range around previous value
        #       maybe: use user-supplied function?
        # todo: create list of calibration experiment expids to submit
        # todo: submit experiments and get rids, then store them in a class dict

        # todo: change status?
        pass

    def _process_calibrations(self):
        """
        todo: document
        """
        # todo: get our class dict which holds calibration results
        # todo: extract optimal value
        # todo: use optimal values to get experiment parameter values

        # todo: check if we have more calibration sets to do, then load them into active calibration queue and submit again
        # todo: otherwise, change status
        pass
