import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: burst counts and burst read


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
    kernel_invariants = {
        "gating_edge", "counting_method"
    }

    def prepare_device(self):
        self.gating_edge = self.get_parameter('gating_edge', group='devices.pmt', override=False)
        # self.time_pmt_gating_mu = self.get_parameter('time_pmt_gating_us', group='devices.pmt', override=False, conversion_function=seconds_to_mu)

        # get default gating edge for counting
        self.counting_method = getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))

        # create counter variable for list timestamping
        self._counter = np.int64(0)

    @kernel(flags={"fast-math"})
    def count(self, time_mu: TInt64) -> TNone:
        """
        Counts the specified gating edges for a given time.
        Arguments:
            time_mu (int64): the time to count for, in machine units (i.e. ns).
        """
        self.counting_method(time_mu)

    @kernel(flags={"fast-math"})
    def timestamp_counts(self, timestamp_mu_list: TArray(TInt64, 1), time_gating_mu: TInt64) -> TNone:
        """
        Record timestamps of incoming counts.
        Does not time out until the passed list (i.e. timestamp_mu_list) is completely filled.
        The passed-in array timestamp_mu_list will be directly modified.

        Arguments:
            timestamp_mu_list   array(int64): an array of timestamps to be filled.
            time_gating_mu      (int64)     : the maximum waiting time (in machine units) for each count.
        """
        # start counting photons
        time_start_mu = now_mu()
        self.input._set_sensitivity(1)

        while self._counter < len(timestamp_mu_list):
            # get count timestamp
            time_photon_mu = self.input.timestamp_mu(now_mu() + time_gating_mu)

            # move timestamped photon into buffer if valid
            if time_photon_mu >= 0:
                timestamp_mu_list[self._counter] = time_photon_mu
                self._counter += 1
            # if count invalid, add slack so we can keep eating counts
            # todo: maybe - should we just terminate if problematic?
            else:
                # todo: maybe make this a timeout slack?
                delay_mu(time_gating_mu)

        # add slack and stop counting
        self.core.break_realtime()
        self.input._set_sensitivity(0)

        # remove timestamping start time
        timestamp_mu_list -= time_start_mu

        # reset loop counter
        # note: we do this here (i.e. at end) to reduce initial overhead where latencies are critical
        self._counter = 0

    @kernel(flags={"fast-math"})
    def clear_inputs(self) -> TNone:
        """
        Clear input count buffers for both the TTLInOut and EdgeCounter objects.
        Note: will consume all slack.
        """
        self.core.break_realtime()  # add slack

        # clear TTL input events
        while self.input.timestamp_mu(now_mu()) != -1:
            delay_mu(5000)
        self.core.break_realtime()

        # clear edge counter input events
        count_events_remaining = self.pmt.fetch_timestamped_count(now_mu())
        while count_events_remaining[0] != -1:
            delay_mu(5000)
            count_events_remaining = self.pmt.fetch_timestamped_count(now_mu())
        self.core.break_realtime()

