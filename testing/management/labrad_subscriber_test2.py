import numpy as np
from artiq.experiment import *

import socket
import labrad
from os import environ
from labrad.wrappers import connectAsync
from twisted.internet.defer import inlineCallbacks, Deferred
from EGGS_labrad.config.multiplexerclient_config import multiplexer_config
TMP_ID = 867234


class LabradSubscriberTest2(EnvExperiment):
    """
    Labrad Subscriber Test2
    Subscribe to intermediate "warning server" on local labrad without
    any personal processing.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("scheduler")

        # set up labrad cxn variables
        # self.LABRADHOST = multiplexer_config.ip
        self.LABRADHOST = environ['LABRADHOST']
        self.LABRADUSERNAME = ""
        self.LABRADPASSWORD = environ['LABRADPASSWORD']

        # create variables for labrad subscription
        self._wm_status_term =  False

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.labrad_prepare()
        self.core.break_realtime()

        self.run_main()
        self.core.break_realtime()
        self.run_cleanup()
        self.core.break_realtime()

    @rpc
    def labrad_prepare(self):
        labrad.thread.startReactor()
        d = Deferred()

        @inlineCallbacks
        def create_connection(msg):
            self.cxn_async = yield connectAsync(
                self.LABRADHOST, port=7682, name="{:s} ({:s})".format("ARTIQ_TEST", socket.gethostname()),
                username=self.LABRADUSERNAME, password=self.LABRADPASSWORD
            )
            self.wm_async = self.cxn_async.warning_server
            self.wm_async.signal__wavemeter_unlock(TMP_ID)
            self.wm_async.addListener(listener=self.update_subscription, source=None, ID=TMP_ID)
            print(msg)

        d.addCallback(create_connection)
        d.callback("\tDEFERRED: FIRED")
        print('\tFINISH LABRAD PREPARE')

    def update_subscription(self, c, signal):
        print("\n\t\tWARNING - CH{:d} UNLOCKED: {:s}".format(*signal))
        self._wm_status_term = True

    @kernel
    def run_main(self):
        self.core.break_realtime()
        while self.check_termination() is False:
            self.core.break_realtime()
            delay_mu(2000000000)

    @rpc
    def run_cleanup(self):
        try:
            self.cxn_async.disconnect()
        except Exception as e:
            print(e)

    @rpc
    def check_termination(self) -> TBool:
        status = self.scheduler.check_termination() or self._wm_status_term
        if self._wm_status_term:
            print("\n\tSTOPPING EXPERIMENT & CANCELLING ALL EXPERIMENTS")
            self.cancel_all_experiments()
        return status

    @rpc
    def cancel_all_experiments(self):
        # get scheduler itinerary
        sched = self.scheduler.get_status()

        # get all experiments in our pipeline
        rid_list = [
            rid
            for rid, exp_dict in sched.items()
            if (rid != self.scheduler.rid) and (exp_dict['pipeline'] == self.scheduler.pipeline_name)
               and (exp_dict['status'] != "running")
        ]
        rid_list.reverse()

        # delete remaining experiments
        for rid in rid_list:
            self.scheduler.delete(rid)
