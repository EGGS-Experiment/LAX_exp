import labrad
import numpy as np
from os import environ
from time import sleep
from artiq.experiment import *

from EGGS_labrad.servers.ARTIQ.artiq_subscriber import ARTIQ_subscriber, ARTIQ_subscriber2
from asyncio import get_event_loop, set_event_loop, Event


class tmptitle(EnvExperiment):
    def build(self):
        self._scheduler = self.get_device("scheduler")
        self._ccb = self.get_device("ccb")
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
        # create necessary values
        self.sbc_freq_list_mhz = np.ones(3) * (self.freq_carrier_mhz - self.freq_secular_khz / 2000.)

        self.scan_range_v = list(self.scan_range_v)
        self.expid_list = list()

        for scan_val_v in self.scan_range_v:
            expid = {
                "file": "experiments_legacy\\parametric_sweep.py",
                "class_name": "ParametricSweep",
                "arguments": {
                    "num_counts": 10000,
                    "ampl_mod_vpp": 0.05,
                    "freq_mod_mhz_list": {
                        "center": 1.206,
                        "span": 0.04,
                        "step": 0.0002,
                        "randomize": True,
                        "seed": None,
                        "ty": "CenterScan"
                    },
                    "dc_micromotion_channel": "V Shim",
                    "dc_micromotion_voltage_v": scan_val_v
                },
                "log_level": 30,
                "repo_rev": "1e672b2317d87dd0e70d49a8f1f1639d3fd217b8"
            }
            self.expid_list.append(expid)


    def run(self):
        # for expid_dj in self.expid_list:
        #     # print(expid_dj["arguments"]["dc_micromotion_voltage_v"])
        #     self._scheduler.submit(pipeline_name="main", expid=expid_dj)
        # connect to master notifications

        try:
            # set up ARTIQ scheduler subscriber
            from asyncio import get_event_loop, set_event_loop, Event
            self.scheduler_subscriber = ARTIQ_subscriber(self._process_scheduler_update, self)
            self.scheduler_subscriber.connect('::1', 3250)

            self.ds_subscriber = ARTIQ_subscriber2(self._process_ds_update, self)
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



        i = 0
        while i < 50:
            sleep(1)
            i+= 1

        self.scheduler_subscriber._subscriber.unsubscribe()
        self.ds_subscriber._subscriber.unsubscribe()


    def _process_scheduler_update(self, mod):
        """
        Checks if any experiments are running and sends a Signal
        to clients accordingly.
        """
        running_exp = None

        # print('\tSCHEDULER UPDATE: {}'.format(mod))

        # check if any experiments are running
        for rid, exp_params in self.struct_holder.backing_store.items():
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
        running_exp = None

        if mod['action'] is not 'init':
            print('\tDS UPDATE: {}'.format(mod))

        # # check if any experiments are running
        # for rid, exp_params in self.struct_holder2.backing_store.items():
        #     run_status = exp_params['status']
        #
        #     # send experiment details to clients if experiment is running
        #     if run_status == 'running':
        #         print('\n\tNEW DS')

        # otherwise, no experiment running, so inform clients
            print('\tIDK DS')

