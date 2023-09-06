import labrad
import numpy as np
from os import environ
from time import sleep

from queue import Queue
from collections import deque

from artiq.experiment import *
from sipyco.sync_struct import Subscriber
from EGGS_labrad.servers.ARTIQ.artiq_subscriber import ARTIQ_subscriber


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

            # set up subscriber objects for the scheduler and the dataset manager
            self.scheduler_subscriber = ARTIQ_subscriber('schedule', self._process_scheduler_update, self)
            self.scheduler_subscriber.connect('::1', 3250)

            self.ds_subscriber = ARTIQ_subscriber('datasets', self._process_ds_update, self)
            self.ds_subscriber.connect('::1', 3250)

            # set up event loop for ARTIQ_subscriber
            loop = get_event_loop()
            stop_event = Event()

            # set up thread to run subscriber event loop in background
            from threading import Thread
            def run_in_background(loop):
                set_event_loop(loop)
                loop.run_until_complete(stop_event.wait())
                loop.close()

            t = Thread(target=run_in_background, args=(loop,))
            t.start()

        except Exception as e:
            print(e)
            print("Unable to connect to ARTIQ Master.")
            raise


        # idk loop
        i = 0
        while i < 50:
            sleep(1)
            i+= 1

        self.scheduler_subscriber._subscriber.unsubscribe()
        self.ds_subscriber._subscriber.unsubscribe()


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


        # print('\tSCHEDULER UPDATE: {}'.format(mod))

        # check if any experiments are running
        for rid, exp_params in self.struct_holder_schedule.backing_store.items():
            run_status = exp_params['status']

            # send experiment details to clients if experiment is running
            if run_status == 'running':
                print('\tNEW EXP RUNNING')

        # otherwise, no experiment running, so inform clients
        # print('\tNO EXP RUNNING')

    def _process_ds_update(self, mod):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        if mod['action'] is not 'init':
            print('\tDS UPDATE: {}'.format(mod))

        # # check if any experiments are running
        # for rid, exp_params in self.struct_holder_datasets.backing_store.items():
        #     run_status = exp_params['status']
        #
        #     # send experiment details to clients if experiment is running
        #     if run_status == 'running':
        #         print('\n\tNEW DS')

        # otherwise, no experiment running, so inform clients
            print('\tIDK DS')


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
