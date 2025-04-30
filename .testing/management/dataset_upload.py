import csv
import labrad
import numpy as np
from os import environ
from artiq.experiment import *


class datasetupload(EnvExperiment):
    """
    dataset upload
    Testing
    """

    def build(self):
        self.setattr_device("core")

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # calib_timestamp = datetime.timestamp(datetime.now())
        # th0 = np.arange(85,137,2)
        # th1 = np.array([0.15625, 0.15625, 0.140625, 0.125, 0.1171875, 0.109375, 0.109375,
        #                 0.109375, 0.1171875, 0.1171875, 0.109375, 0.109375, 0.109375,
        #                 0.1171875, 0.125, 0.125, 0.125, 0.1328125, 0.140625, 0.140625,
        #                 0.15625, 0.171875, 0.203125, 0.25, 0.28125, 0.34375]) * 100

        # # print(np.array([th0,th1]))
        # self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', np.array([th0, th1]).transpose(), broadcast=True, persist=True)
        # self.set_dataset('calibration.temperature.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

    def prepare(self):
        fname= 'C:\\Users\\EGGS1\\Documents\\res1.csv'
        with open(fname, newline='') as f:
            reader = csv.reader(f)
            th0 = []
            for row in reader:
                th0.append([float(val) for val in row])

        th0=np.array(th0)
        self.set_dataset('calibration.photodiode.responsivity.DET36A2', th0, broadcast=True,persist=True)

    def run(self):
        pass

    def analyze(self):
        pass
