import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class RPCTimeProfiler(EnvExperiment):
    """
    RPC Time Profiler
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device('urukul2_ch1')

    def prepare(self):
        self.repetitions = 1000
        self.time_delay_mu = self.core.seconds_to_mu(3 * ms)

        self.set_dataset("results", np.zeros(self.repetitions, dtype=np.int64))
        self.setattr_dataset("results")
        self.idx_results = 0

    @rpc(flags={"async"})
    def update_results(self, time_mu: TInt64):
        self.mutate_dataset('results', self.idx_results, time_mu)
        self.idx_results += 1

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()
        self.urukul2_ch1.set(110*MHz, amplitude=0.5)
        self.core.break_realtime()

        for i in range(self.repetitions):
            self.core.break_realtime()

            # simulate pulse sequence
            self.urukul2_ch1.sw.on()
            delay_mu(self.time_delay_mu)
            self.urukul2_ch1.sw.off()

            # RPC call
            time_start_mu = self.core.get_rtio_counter_mu()
            self.rpc_call()
            time_stop_mu = self.core.get_rtio_counter_mu()
            self.update_results(time_stop_mu - time_start_mu)

        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @rpc
    def rpc_call(self) -> TNone:
        return

