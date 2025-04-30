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

        self.setattr_argument("psd_units",          BooleanValue(default=True))
        self.setattr_argument("span_freq_khz_list", Scannable(
                                                            default=[
                                                                ExplicitScan([12.8, 103.]),
                                                            ],
                                                            global_min=0.01, global_max=105, global_step=1,
                                                            unit="kHz", scale=1, precision=3
                                                        ))

    def prepare(self):
        """
        todo: document
        """
        # ensure display num is either 0 or 1
        self.display_num = round(self.display_num)

        # save relevant parameters as dataset (for convenience later on)
        self.set_dataset('config', self.config_string, broadcast=False)
        self.set_dataset('display_num', self.display_num, broadcast=False)
        self.set_dataset('calibration_factor', self.calibration_factor, broadcast=False)
        self.set_dataset('psd_units', self.calibration_factor, broadcast=False)

        # counters
        self._idx_results = 0   # result dataset counter
        self._idx_failure = 0   # failure counter

        # magic numbers
        self.max_failures = 5       # max number of failures
        self._poll_delay_s = 2.5    # poll delay document


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

        # store/set up important parameters
        # set up PSD units
        self.dev.write("PSDU {:d},{:d}".format(self.display_num, self.psd_units))
        # todo: store input config: grounding, coupling, range
        # todo: store measurement: type, units
        # todo: store averaging: type, num_avg

    def get_trace(self, freq_span_khz: TFloat) -> TNone:
        """
        Get 2D trace values.
        Arguments:
            freq_span_khz: the frequency span to use for the trace (in kHz).
        """
        # set recording span & start data acquisition
        self.dev.write("*CLS; FSPN {:d},{:f}; STRT".format(2, freq_span_khz * kHz))

        # poll until avgs completed
        completion_status = bool(int(self.dev.query("DSPS? {:d}".format(round(1 + self.display_num * 8)))))
        while completion_status is False:
            sleep(self._poll_delay_s)

            try:
                completion_status = bool(int(self.dev.query("DSPS? {:d}".format(round(1 + self.display_num * 8)))))

            except Exception as e:
                print("Failure #{:d}: {}".format(self._idx_failure, e))
                self._idx_failure += 1
                if self._idx_failure > self.max_failures:
                    raise ValueError("Maximum number of failures reached.")

        # get trace data
        try:
            # get frequency start/stop to create x-axis array in Hz
            freq_start_hz = float(self.dev.query("FSTR? {:d}".format(self.display_num)).strip())
            freq_stop_hz = float(self.dev.query("FEND? {:d}".format(self.display_num)).strip())
            num_freq_points = int(self.dev.query("DSPN? {:d}".format(self.display_num)).strip())
            xvals_freq_hz = np.linspace(freq_start_hz, freq_stop_hz, num_freq_points)

            # get trace values and convert to float (units should be dBm)
            yvals_ampl_dbm = self.dev.query("DSPY? {:d}".format(self.display_num))
            yvals_ampl_dbm = np.array([float(str_val) for str_val in yvals_ampl_dbm.split(",")])
            # update results w/calibration factor
            yvals_ampl_dbm += self.calibration_factor

        except Exception as e:
            self._idx_failure += 1
            print("Failure #{:d}: {}".format(self._idx_failure, e))
            if self._idx_failure > self.max_failures:
                raise ValueError("Maximum number of failures reached.")

        else:
            # save results in 2D dataset
            self.set_dataset(
                "results_{:d}".format(self._idx_results),
                np.transpose([xvals_freq_hz, yvals_ampl_dbm]),
                broadcast=False
            )
            self._idx_results += 1

    def cleanup_device(self) -> TNone:
        """
        Clean up device settings after we're done.
        """
        # set PSD units back to normal
        self.dev.write("PSDU {:d},{:d}".format(self.display_num, False))


    def run(self) -> TNone:
        """
        Main sequence.
        """
        try:
            # connect to device & set up
            self.prepare_pyvisa()

            # get trace
            for freq_span_khz in self.span_freq_khz_list:
                self.get_trace(freq_span_khz)

            # clean up
            self.cleanup_device()

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
