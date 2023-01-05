from artiq.experiment import *

class ScheduleTest(EnvExperiment):
    def build(self):
        self._scheduler = self.get_device("scheduler")

    def run(self):
        new_expid = {
            "file": "repository/tmp2.py",
            "class_name": "testarg12",
            "arguments": {},
            "log_level": 30,
        }

        self._scheduler.submit(pipeline_name="main", expid=new_expid)
