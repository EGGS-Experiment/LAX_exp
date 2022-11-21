import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class Testing(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl_counter3")

        self.time0 = self.core.seconds_to_mu(1 * ms)
        self.time1 = self.core.seconds_to_mu(10 * us)
        self.timed = self.core.seconds_to_mu(5 * us)

        self.set_dataset("tmp", [])
        self.setattr_dataset("tmp")

        # self.set_dataset('dds_tickle_channel', 3, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    @kernel
    def run(self):
        self.core.reset()

        self.ttl0.input()
        #self.ttl3.input()
        self.core.break_realtime()

        self.ttl_counter3.gate_rising_mu(self.time0)
        print(self.ttl_counter3.fetch_count())
        #self.append_to_dataset('tmp', self.ttl_counter3.fetch_count())

        # for i in range(100):
        #
        #     # th1
        #     # self.ttl0._set_sensitivity(1)
        #     # delay_mu(self.time0)
        #     # time_input_mu = self.ttl0.timestamp_mu()
        #
        #     # th0
        #     time_end_mu = self.ttl0.gate_rising_mu(self.time0)
        #     time_input_mu = self.ttl0.timestamp_mu(time_end_mu)
        #
        #     # th2
        #     if time_input_mu > 0:
        #         at_mu(time_input_mu)
        #         delay_mu(self.timed)
        #         time_end_mu2 = self.ttl3.gate_rising_mu(self.time1)
        #         time_input_mu2 = self.ttl3.timestamp_mu(time_end_mu2)
        #         #self.ttl3.count(time_end_mu2)
        #
        #         #self.ttl0.count(time_end_mu)
        #         self.core.break_realtime()
        #
        #         self.append_to_dataset("tmp", [time_input_mu, time_input_mu2])
        #         self.core.break_realtime()

    def analyze(self):
        # self.tmp = np.array([self.core.mu_to_seconds(val[1] - val[0]) for val in self.tmp])
        # print(self.tmp)
        # print(len(self.tmp))
        self.tmp = np.array(self.tmp)
        print(self.tmp)
