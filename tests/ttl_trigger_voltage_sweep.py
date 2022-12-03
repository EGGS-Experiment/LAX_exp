import labrad
import numpy as np

from os import environ
from artiq.experiment import *


class TTLTriggerVoltageSweep(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0")
        self.setattr_device("ttl3")

        self.repetitions = 20000

        self.time0 = self.core.seconds_to_mu(1 * ms)
        self.time1 = self.core.seconds_to_mu(10 * us)
        self.timed = self.core.seconds_to_mu(5 * us)

        self.set_dataset("ttl_trigger", [])
        self.setattr_dataset("ttl_trigger")

        self.set_dataset("ttl_trigger_processed", np.zeros(self.repetitions))
        self.setattr_dataset("ttl_trigger_processed")

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc = self.cxn.dc_server


    def prepare(self):
        self.voltage_channel = 20
        self.voltage_list_v = np.arange(40, 80 + 1, 1)
        self.setattr_dataset('xArr', self.voltage_list_v)


    @kernel
    def run(self):
        self.core.reset()

        # set ttl direction
        self.ttl0.input()
        self.ttl3.input()
        self.core.break_realtime()


        # MAIN LOOP
        # sweep voltage
        for voltage_val in self.voltage_list_v:

            # set voltage
            self.voltage_set(self.voltage_channel, voltage_val)
            self.core.break_realtime()

            # get photon counts
            for i in range(self.repetitions):
                self.core.break_realtime()

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
                    self.update_dataset(voltage_val, time_input_mu, time_input_mu2)
                    self.core.break_realtime()


    @rpc(flags={"async"})
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        print('voltage set: {}'.format(voltage_set_v))


    @rpc(flags={"async"})
    def update_dataset(self, voltage_v, time_start_mu, time_stop_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('ttl_trigger', [voltage_v, time_start_mu, time_stop_mu])


    def analyze(self):
        pass
        # get photon correlation time
        # self.ttl_trigger_processed = np.array([self.core.mu_to_seconds(val[1] - val[0]) for val in self.ttl_trigger])
        # remove constant offset
        #self.ttl_trigger_processed -= np.amin(self.ttl_trigger_processed)
        #print(self.ttl_trigger_processed)
