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
        for expid_dj in self.expid_list:
            # print(expid_dj["arguments"]["dc_micromotion_voltage_v"])
            self._scheduler.submit(pipeline_name="main", expid=expid_dj)
