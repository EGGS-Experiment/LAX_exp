import labrad
import numpy as np
from os import environ

from artiq.experiment import *


class ScheduleTest(EnvExperiment):
    def build(self):
        self._scheduler = self.get_device("scheduler")
        self.setattr_argument("scan_range_v",                           Scannable(
                                                                            default=CenterScan(40, 80, 1, randomize=True),
                                                                            global_min=0, global_max=300, global_step=1,
                                                                            unit="V", scale=1, ndecimals=3
                                                                        ))


    def prepare(self):
        # self.scan_range_v = list(self.scan_range_v)
        self.scan_range_v = [1,2]
        self.expid_list = list()

        for scan_val_v in self.scan_range_v:
            expid = {
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
            }
            self.expid_list.append(expid)


    def run(self):
        for expid_dj in self.expid_list:
            # print(expid_dj["arguments"]["dc_micromotion_voltage_v"])
            self._scheduler.submit(pipeline_name="main", expid=expid_dj)
