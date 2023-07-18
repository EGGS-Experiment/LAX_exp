from numpy import array
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


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

    PMT_KEY = "pmt_counts"


    def prepare_device(self):
        self.gating_edge =                              self.get_parameter('gating_edge', group='pmt', override=False)
        # self.time_pmt_gating_mu =                       self.get_parameter('time_pmt_gating_us', group='pmt', override=False, conversion_function=seconds_to_mu)

        # get default gating edge for counting
        self.counting_method =                          getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.set_dataset(self.PMT_KEY, [])
        self.kernel_invariants.add('counting_method')

    @kernel(flags={"fast-math"})
    def count(self, time_mu: TInt64):
        """
        Counts the specified gating edges for a given time.
        Arguments:
            time_mu (int64): the time to count for, in machine units (i.e. ns).
        """
        self.counting_method(time_mu)

    @kernel(flags={"fast-math"})
    def timestamp_counts(self, num_counts: TInt32, time_gating_mu: TInt64) -> TArray(TInt64, 1):
        """
        Timestamp the incoming counts.
        Does not time out until we read finish reading the given number of counts.
        Arguments:
            num_counts      (int)   : the number of counts to read.
            time_gating_mu  (int64) : the maximum waiting time (in machine units) for each count.
        Returns:
                            list(int64): the list of count timestamps (in machine units).
        """
        # set counting variables
        counter = 0
        timestamp_mu_list = [0] * num_counts
        self.core.break_realtime()

        # start counting photons
        self.input._set_sensitivity(1)
        while counter < num_counts:

            # get count timestamp
            time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

            # move timestamped photon into buffer if valid
            if time_photon_mu >= 0:
                timestamp_mu_list[counter] = time_photon_mu
                counter += 1

        # stop counting
        self.core.break_realtime()
        self.input._set_sensitivity(0)

        # return timestamps
        return array(timestamp_mu_list)

    # @kernel(flags={"fast-math"})
    # def get_timestamp(self) -> TInt64:
    #     """
    #     # tmp rename
    #     Timestamp the incoming counts.
    #     Does not time out until we read finish reading the given number of counts.
    #     Arguments:
    #         num_counts      (int)   : the number of counts to read.
    #         time_gating_mu  (int64) : the maximum waiting time (in machine units) for each count.
    #     Returns:
    #                         list(int64): the list of count timestamps (in machine units).
    #     """
    #     # set counting variables
    #     counter = 0
    #     timestamp_mu_list = [0] * num_counts
    #     # self.core.break_realtime()
    #
    #     # start counting photons
    #     self.input._set_sensitivity(1)
    #     while counter < num_counts:
    #
    #         # get count timestamp
    #         time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)
    #
    #         # move timestamped photon into buffer if valid
    #         if time_photon_mu >= 0:
    #             timestamp_mu_list[counter] = time_photon_mu
    #             counter += 1
    #
    #     # stop counting
    #     self.core.break_realtime()
    #     self.input._set_sensitivity(0)
    #
    #     # return timestamps
    #     return array(timestamp_mu_list)

    @kernel(flags={"fast-math"})
    def timestamp_counts_fixed_time(self, time_interval_mu: TInt32, time_gating_mu: TInt64) -> TArray(TInt64, 1):
        """
        Timestamp the incoming counts.
        Does not time out until we read finish reading the given number of counts.
        Arguments:
            time_interval_mu      (int)   : machine units to read processes
            time_gating_mu  (int64) : the maximum waiting time (in machine units) for each count.
        Returns:
                            list(int64): the list of count timestamps (in machine units).
        """
        # start counting photons
        self.input.gate_rising_mu(time_interval_mu)

        # get count timestamp
        time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

        # move timestamped photon into buffer if valid
        if time_photon_mu >= 0:
            while time_photon_mu >= 0:
                self.append_to_dataset(self.PMT_KEY, time_photon_mu)
                time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

        # stop counting
        self.core.break_realtime()
        self.input._set_sensitivity(0)
        timestamp_mu_list = self.get_dataset(self.PMT_KEY)
        # return timestamps
        return array(timestamp_mu_list)
