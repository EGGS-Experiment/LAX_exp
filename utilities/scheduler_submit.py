from artiq.experiment import *

class SchedulerExperimentScan(EnvExperiment):
    def build(self):
        self._scheduler = self.get_device("scheduler")
        # todo: add exp name/file/whatever as arg
        # todo: set argument - dict of standard args
        # todo: set argument - argument to be scanned

    def prepare(self):
        pass
        # todo: create array of expids

    def run(self):
        # todo: wait until one experiment has finished before submitting the next one
        new_expid = {
            "file": "LAX_exp/tmp2.py",
            "class_name": "testarg12",
            "arguments": {},
            "log_level": 30,
        }

        self._scheduler.submit(pipeline_name="main", expid=new_expid)
