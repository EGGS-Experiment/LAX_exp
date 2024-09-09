import numpy as np
from artiq.experiment import *

from time import sleep


class NetanalMeasure(EnvExperiment):
    """
    Utility: Netanal Measure

    Store a trace from the network analyzer.
    """

    def build(self):
        pass

    def prepare(self):
        _CONFIG = "peterbox_s11"

        # connect to labrad
        import labrad
        self.cxn=labrad.connect()

        # set up network analyzer device
        self.na = self.cxn.network_analyzer_server
        self.na.select_device()
        self.num_points = self.na.sweep_points(1601)
        self.freq_start = self.na.frequency_start()
        self.freq_stop = self.na.frequency_stop()

        # autoscale y-axes
        self.na.gpib_write('DISP:WIND1:TRAC:Y:AUTO ONCE')
        self.na.gpib_write('DISP:WIND2:TRAC:Y:AUTO ONCE')

        # set up data structures and datasets
        self._results0 = np.linspace(self.freq_start, self.freq_stop, self.num_points)
        self._results1 = np.zeros(self.num_points)

        self.set_dataset('config', _CONFIG, broadcast=False)
        self.set_dataset('freq_start_hz', self.freq_start, broadcast=False)
        self.set_dataset('freq_stop_hz', self.freq_stop, broadcast=False)
        self.set_dataset('num_points', self.num_points, broadcast=False)

    def run(self):
        try:
            # lock device to prevent communication interruptions
            self.na.lock_device()
            res10 = self.na.gpib_query('TRAC? CH1FDATA')
            res11 = np.array([float(strval) for strval in res10.split(',')])
            self._results1 = res11
            # print('\trun 1 success')
            sleep(2.)

            # store channel 2
            # print('prepare to run (2)')
            # sleep(2.)
            res20 = self.na.gpib_query('TRAC? CH2FDATA')
            res21 = np.array([float(strval) for strval in res20.split(',')])
            self._results2 = res21
            # print('\trun 2 success')
            # sleep(2.)
        except Exception as e:
            print("Error: unable to save network analyzer traces.")
            print(repr(e))
        finally:
            self.na.release_device()
            self.cxn.disconnect()

    def analyze(self):
        # save trace from channel 1
        res_fin1 = np.array([self._results0, self._results1]).transpose()
        self.set_dataset('results1', res_fin1)
        # print(res_fin1)

        # save trace from channel 2
        res_fin2 = np.array([self._results0, self._results2]).transpose()
        self.set_dataset('results2', res_fin2)
        # print(res_fin2)
