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
    devices ={
        'input':    'ttl0'
    }


    def prepare_device(self):
        self.gating_edge =                              self.get_parameter('gating_edge', group='pmt', override=False)
        # self.time_pmt_gating_mu =                       self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu)

        # get default gating edge for counting
        self.counting_method =                          getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')

        # todo: create counter variable for timestamping

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
        # todo: write that it modifies the list object passed to it

        Arguments:
            # todo: modify and fix
            num_counts      (int)   : the number of counts to read.
            time_gating_mu  (int64) : the maximum waiting time (in machine units) for each count.
        Returns:
                            list(int64): the list of count timestamps (in machine units).
        """
        # set counting variables
        counter = 0

        # start counting photons
        self.input._set_sensitivity(1)
        while counter < len(timestamp_mu_list):

            # get count timestamp
            time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

            # move timestamped photon into buffer if valid
            if time_photon_mu >= 0:
                timestamp_mu_list[counter] = time_photon_mu
                counter += 1

        # stop counting
        with parallel:
            self.core.break_realtime()
            self.input._set_sensitivity(0)
