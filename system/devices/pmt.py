import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: aggressive optimization


class PMTCounter(LAXDevice):
    """
    Device: PMT photon counter

    Wrapper for the TTL edge_counter and ttl_input objects that read in PMT counts.
    """
    name = "pmt"
    core_device = ('pmt', 'ttl0_counter')
    devices = {
        'input':    'ttl0'
    }


    def prepare_device(self):
        self.gating_edge =                              self.get_parameter('gating_edge', group='pmt', override=False)
        # self.time_pmt_gating_mu =                       self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu)

        # get default gating edge for counting
        self.counting_method =                          getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')

        # create counter variable for list timestamping
        self.counter = np.int64(0)

    @kernel(flags={"fast-math"})
    def count(self, time_mu: TInt64):
        """
        Counts the specified gating edges for a given time.
        Arguments:
            time_mu (int64): the time to count for, in machine units (i.e. ns).
        """
        self.counting_method(time_mu)

    @kernel(flags={"fast-math"})
    def timestamp_counts(self, timestamp_mu_list: TArray(TInt64, 1), time_gating_mu: TInt64):
        """
        Timestamp the incoming counts.
        Does not time out until we read finish reading the given number of counts.
        The passed-in array timestamp_mu_list will be directly modified.

        Arguments:
            timestamp_mu_list   array(int64): an array of timestamps to be filled.
            time_gating_mu      (int64)     : the maximum waiting time (in machine units) for each count.
        """
        # start counting photons
        time_start_mu = now_mu()
        self.input._set_sensitivity(1)

        while self.counter < len(timestamp_mu_list):
            # get count timestamp
            time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

            # move timestamped photon into buffer if valid
            if time_photon_mu >= 0:
                timestamp_mu_list[self.counter] = time_photon_mu
                self.counter += 1

        # add slack and stop counting
        self.core.break_realtime()
        self.input._set_sensitivity(0)

        # remove timestamping start time
        timestamp_mu_list -= time_start_mu

        # reset loop counter
        # note: we do this here (i.e. at end) to reduce initial overhead where latencies are critical
        self.counter = 0
