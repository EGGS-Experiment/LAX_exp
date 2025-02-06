import numpy as np
from artiq.experiment import *

import pyvisa
from time import sleep


class FFTAnalyzerTrace(EnvExperiment):
    """
    Tool: FFT Analyzer Trace

    Store a trace from the SR7xx FFT analyzers.
    """

    def build(self):
        self.setattr_argument("device_key",         StringValue(default="GPIB1::10"))
        self.setattr_argument("display_num",        NumberValue(default=0, precision=0, step=1, min=0, max=1))
        self.setattr_argument("config_string",      StringValue(default="N/A"))
        self.setattr_argument("calibration_factor", NumberValue(default=5., precision=2, step=10, min=-100, max=100))

    def prepare(self):
        """
        todo: document
        """
        # ensure display num is either 0 or 1
        self.display_num = round(self.display_num)

        # set relevant parameters as dataset
        self.set_dataset('config', self.config_string, broadcast=False)
        self.set_dataset('display_num', self.display_num, broadcast=False)
        self.set_dataset('calibration_factor', self.calibration_factor, broadcast=False)


    """
    PYVISA METHODS
    """
    def prepare_pyvisa(self) -> TNone:
        """
        Open pyvisa, connect to device, and configure communication.
        """
        # instantiate objects and get resouce list
        self.rm = pyvisa.ResourceManager()
        lr = self.rm.list_resources()

        # look for target device
        try:
            dev_name = next(dev_name for dev_name in lr if self.device_key in dev_name)
        except StopIteration:
            raise ValueError("Unable to find target GPIB device in pyvisa.")

        # open target device and configure
        self.dev = self.rm.open_resource(dev_name)
        self.dev.write_termination = '\n'   # default is '\r\n', which is invalid for SR7xx
        self.dev.timeout = 10000            # in ms, default is 3s

        # todo: store important parameters
        # todo: store PSD units?
        # todo: store linear vs vector averaging
        # todo: store units
        # todo: store measurement type

    def get_trace(self) -> TNone:
        """
        Get 2D trace values
        """
        # get frequency start/stop
        self.freq_start_hz =    float(self.dev.query("FSTR? {:d}".format(self.display_num)).strip())
        self.freq_stop_hz =     float(self.dev.query("FEND? {:d}".format(self.display_num)).strip())
        self.num_freq_points =  int(self.dev.query("DSPN? {:d}".format(self.display_num)).strip())
        # create x-axis array in Hz
        xvals_freq_hz = np.linspace(self.freq_start_hz, self.freq_stop_hz, self.num_freq_points)

        # get trace values and convert to float (units should be dBm)
        yvals_ampl_dbm = self.dev.query("DSPY? {:d}".format(self.display_num))
        yvals_ampl_dbm = np.array([float(str_val) for str_val in yvals_ampl_dbm.split(",")])
        # update results w/calibration factor
        yvals_ampl_dbm += self.calibration_factor

        # save results in 2D dataset
        self.set_dataset(
            "results",
            np.transpose([xvals_freq_hz, yvals_ampl_dbm]),
            broadcast=False
        )

    def run(self) -> TNone:
        """
        todo: document
        """
        try:
            # connect to device & set up
            self.prepare_pyvisa()

            # get trace
            self.get_trace()

        except Exception as e:
            # print error message
            print("Error: {}".format(e))
            print(repr(e))

        finally:
            # make sure to close device!
            if hasattr(self, 'dev'):
                self.dev.close()

    def analyze(self):
        pass
