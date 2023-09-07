import labrad
import numpy as np
from os import environ

from artiq.experiment import *

from time import sleep
from queue import Queue
from collections import deque

# import asyncio
from sipyco.sync_struct import Subscriber
from EGGS_labrad.servers.ARTIQ.artiq_subscriber import ARTIQ_subscriber
from artiq.frontend.artiq_client import _run_subscriber


class tmptitle(EnvExperiment):
    def build(self):
        self.ccb = self.get_device("ccb")
        self.scheduler = self.get_device("scheduler")
        self.setattr_argument("scan_range_v",                           Scannable(
                                                                            default=CenterScan(40, 80, 1, randomize=True),
                                                                            global_min=0, global_max=300, global_step=1,
                                                                            unit="V", scale=1, ndecimals=3
                                                                        ))

        # base params - laser scan
        self.freq_carrier_mhz =     103.209
        self.freq_secular_khz =     1088.0
        self.freq_ls_scan_khz =     25.0

        # base params - sbc
        self.sbc_freq_span_khz =    4.
        self.sbc_freq_points =      5.


    def prepare(self):
        # create necessary data structures
        self.experiment_queue =         deque()
        self.submitted_experiments =    deque()
        self.calibration_results =      dict()

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
            from asyncio import get_event_loop, set_event_loop, Event
            from sipyco.asyncio_tools import atexit_register_coroutine
            #
            # # set up subscriber objects for the scheduler and the dataset manager
            # self.scheduler_subscriber = ARTIQ_subscriber('schedule', self._process_scheduler_update, self)
            # self.scheduler_subscriber.connect('::1', 3250)
            #
            # self.ds_subscriber = ARTIQ_subscriber('datasets', self._process_ds_update, self)
            # self.ds_subscriber.connect('::1', 3250)
            #
            # # set up event loop for ARTIQ_subscriber
            # loop = get_event_loop()
            # stop_event = Event()
            #
            # # set up thread to run subscriber event loop in background
            # from threading import Thread
            # def run_in_background(loop):
            #     set_event_loop(loop)
            #     loop.run_until_complete(stop_event.wait())
            #     loop.close()
            #
            # t = Thread(target=run_in_background, args=(loop,))
            # t.start()
            _scheduler_struct = dict()
            def _update_scheduler(x):
                print('\t\t\t\t\t\t_update_scheduler_called')
                _scheduler_struct.clear()
                _scheduler_struct.update(x)
                return _scheduler_struct

            _dataset_struct = dict()
            def _update_datasets(x):
                print('\t\t\t\t\t\t_update_datasets called')
                _dataset_struct.clear()
                _dataset_struct.update(x)
                return _dataset_struct


            # set up event loop for subscribers
            loop = get_event_loop()
            stop_event = Event()

            # create subscribers
            # self.scheduler_subscriber = Subscriber('schedule', _update_scheduler, self._process_scheduler_update)
            # self.dataset_subscriber = Subscriber('datasets', _update_datasets, self._process_datasets_update)
            self.scheduler_subscriber = Subscriber('schedule', _update_scheduler, lambda mod: self._process_scheduler_update(_scheduler_struct))
            self.dataset_subscriber = Subscriber('datasets', _update_datasets, lambda mod: self._process_dataset_update(_dataset_struct))
            # connect subscribers
            loop.run_until_complete(self.scheduler_subscriber.connect('::1', 3250))
            loop.run_until_complete(self.dataset_subscriber.connect('::1', 3250))
            # ensure subscribers close connection upon exit
            atexit_register_coroutine(self.scheduler_subscriber.close)
            atexit_register_coroutine(self.dataset_subscriber.close)

            # set up thread to run subscriber event loop in background
            # from threading import Thread
            # def run_in_background(loop):
            #     # set thread event loop as given
            #     set_event_loop(loop)
            #     loop.run_until_complete(stop_event.wait())
            #     loop.close()
            #
            # t = Thread(target=run_in_background, args=(loop,))
            # t.start()

            # ruh event loop indefinitely
            loop.run_until_complete(stop_event.wait())
            loop.close()

        except Exception as e:
            print('\n\t\t error: {}'.format(e))
            raise


        # idk loop
        i = 0
        while i < 1000:
            sleep(5)
            i+= 1
            print('h_loop_{}'.format(i))
            print('\t\tsched: {}'.format(list(_scheduler_struct.keys())))
            print('\t\tkeys: {}'.format(list(_dataset_struct.keys())[-1]))


    '''
    Subscriber idk todo document
    '''

    def _process_scheduler_update(self, mod):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        # todo: update holding structure that contains RIDs
        # todo: check if

        print('\tSCHEDULER UPDATE: {}'.format(mod))
        # check if any experiments are running
        # for rid, exp_params in self._scheduler_struct.backing_store.items():
        #     run_status = exp_params['status']
        #
        #     # send experiment details to clients if experiment is running
        #     if run_status == 'running':
        #         print('\tNEW EXP RUNNING')

    def _process_datasets_update(self, mod):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        # print('\tDS UPDATE: {}'.format(mod))
        try:
            pass
            if mod['action'] != 'init':
                if (mod['key'] != 'progress'):
                    pass
                    if (mod['key'] != 'management.completion_pct'):
                        if int(mod['value'][1] * 10) % 5:
                            print('\tDS UPDATE: {}'.format(mod))
                            print('\t\tkeys: {}'.format(list(self._dataset_struct.keys())[-1]))
        except Exception as e:
            print('\n\tds err: {}'.format(e))
        # if mod['action'] is not 'init':
        #     print('\tDS UPDATE: {}'.format(mod))

        # # check if any experiments are running
        # for rid, exp_params in self.struct_holder_datasets.backing_store.items():
        #     run_status = exp_params['status']
        #
        #     # send experiment details to clients if experiment is running
        #     if run_status == 'running':
        #         print('\n\tNEW DS')

        # otherwise, no experiment running, so inform clients
        # print('\tIDK DS')


    '''
    idk logic i guess
    '''

    def _submit_calibrations(self):
        # todo:
        pass

    def _process_calibrations(self):
        pass

    def _submit_experiments(self, num_exps, param_dict):
        # todo: loop over num_exps
        # todo: get exp from tosubmit experiment expid storage
        # todo: updaet expid dict with our new calibrated parameters
        # todo: submit updated expid to scheduler
        # todo: add RID of submitted exp to rid storage to allow scheduler processor to check completion


        # for expid_dj in self.expid_list:
        #     # print(expid_dj["arguments"]["dc_micromotion_voltage_v"])
        #     self.scheduler.submit(pipeline_name="main", expid=expid_dj)
        pass
