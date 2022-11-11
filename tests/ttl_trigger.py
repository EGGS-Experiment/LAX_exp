import numpy as np
from artiq.experiment import *


class TTLTrigger(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl3")

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

        # set ttl direction
        self.ttl0.input()
        self.ttl3.input()
        self.core.break_realtime()


        # MAIN LOOP
        for i in range(100):

            # wait for event
            time_end_mu = self.ttl0.gate_rising_mu(self.time0)
            time_input_mu = self.ttl0.timestamp_mu(time_end_mu)

            # check if event has fired
            if time_input_mu > 0:

                # set RTIO time and add slack
                at_mu(time_input_mu)
                delay_mu(self.timed)

                # get timestamp of RF event
                time_end_mu2 = self.ttl3.gate_rising_mu(self.time1)
                time_input_mu2 = self.ttl3.timestamp_mu(time_end_mu2)

                # close input gating
                self.ttl3.count(time_end_mu2)
                self.ttl0.count(time_end_mu)
                self.core.break_realtime()

                # add data to dataset
                self.append_to_dataset("tmp", [time_input_mu, time_input_mu2])
                self.core.break_realtime()

    def analyze(self):
        # get photon correlation time
        self.tmp = np.array([self.core.mu_to_seconds(val[1] - val[0]) for val in self.tmp])
        # remove constant offset
        self.tmp -= np.amin(self.tmp)
        print(self.tmp)
