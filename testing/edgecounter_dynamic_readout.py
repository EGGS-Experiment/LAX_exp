import math
import numpy as np
from artiq.experiment import *
from artiq.coredevice.rtio import rtio_output, rtio_input_timestamp, rtio_input_data
from artiq.coredevice.edge_counter import (CONFIG_COUNT_RISING, CONFIG_COUNT_FALLING,
                                           CONFIG_SEND_COUNT_EVENT, CONFIG_RESET_TO_ZERO)


class TTLDynamicTest(EnvExperiment):
    """
    TTLDynamicTest
    idk 2 TTLDynamicTest
    """
    kernel_invariants = {
        # devices
        "ttl", "ttl_chan",

        # timing
        "time_readout_us", "time_readout_mu", "time_bin_us", "time_bin_mu",

        # config/binning
        "repetitions", "num_sub_bins",

        # count rates & probabilities
        "count_rate_bright", "count_rate_dark", "count_rate_bright_bin", "count_rate_dark_bin", "max_counts_bin",
        "likelihood_bright_n", "likelihood_dark_n", "error_threshold", "error_threshold_frac",
        "_dist_sigma_max"
    }

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl0_counter")

        self.repetitions =      1000
        self.time_readout_us =  3000
        self.time_bin_us =      500

        self.count_rate_bright =    120     # per 3ms
        self.count_rate_dark =      30      # per 3ms
        self.error_threshold =      0.01    # total error fraction

        self._dist_sigma_max =  10  # \sigma above mean for max counts

    def prepare(self):
        self.ttl =      self.get_device("ttl0_counter")
        self.ttl_chan = self.ttl.channel << 8

        self.time_bin_mu =      self.core.seconds_to_mu(self.time_bin_us * us)
        self.num_sub_bins =     self.time_readout_us // self.time_bin_us
        self.time_readout_mu =  self.core.seconds_to_mu((self.time_bin_us * us) * self.num_sub_bins)

        self.count_rate_bright_bin =    np.int32(self.count_rate_bright * ((self.time_bin_us * us) / (3. * ms)))
        self.count_rate_dark_bin =      np.int32(self.count_rate_dark * ((self.time_bin_us * us) / (3. * ms)))

        # reprocess error threshold for faster calculation (avoids division
        self.error_threshold_frac =     self.error_threshold / (1. - self.error_threshold)


        # set max counts as 10\sigma above distribution parameter
        self.max_counts_bin = round(
            (self.count_rate_bright_bin + self.count_rate_dark_bin) +
            self._dist_sigma_max * math.sqrt(self.count_rate_bright_bin + self.count_rate_dark_bin)
        )

        # precalculate factorials
        # todo: use normal factorial for n<=10, and use stirling's approx for n > 10
        factorial_arr = np.array([math.factorial(n)
                                  for n in range(0, self.max_counts_bin + 1)])

        # precalculate all relevant likelihoods
        self.likelihood_bright_n =  (
                math.exp(-1. * self.count_rate_bright_bin) *
                (self.count_rate_bright_bin ** np.arange(0, self.max_counts_bin + 1)) /
                factorial_arr
        )
        self.likelihood_dark_n =    (
                math.exp(-1. * self.count_rate_dark_bin) *
                (self.count_rate_dark_bin ** np.arange(0, self.max_counts_bin + 1)) /
                factorial_arr
        )

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset
        self.core.reset()
        delay_mu(1000000)   # 1ms

        # todo: add error handling

        for i in range(self.repetitions):
            # simulate pulse sequence
            self.core.break_realtime()
            delay_mu(1000000)   # 1ms

            # dynamic readout
            self._tmpfunc()

            # add slack/finish up
            self.core.break_realtime()
            delay_mu(125000)    # 125us

    @kernel(flags={"fast-math"})
    def _tmpfunc(self) -> TTuple([TInt32, TInt32, TInt32]):
        """
        todo: document
        Returns: a tuple of (ion_state, total_counts, elapsed_bins).
        """
        # get reference time
        time_start_mu = now_mu()

        # instantiate variables
        ion_state =     -1
        total_counts =  0
        bin_counter =   0
        p_b =           0   # likelihood bright
        p_d =           0   # likelihood dark

        # start initial sub-bin
        at_mu(time_start_mu)
        rtio_output(self.ttl_chan, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)

        # dynamically process each sub-bin
        while bin_counter < self.num_sub_bins:
            # # close previous sub-bin
            at_mu(now_mu() + self.time_bin_mu)
            rtio_output(self.ttl_chan, CONFIG_SEND_COUNT_EVENT)
            # # prepare next sub-bin
            delay_mu(8)
            rtio_output(self.ttl_chan, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)

            # todo: can i send count event while still gating/eating counts?
            # todo: if I send count event and reset_to_zero, does that clear the counts from the count event?

            # eat counts and update loop
            counts_tmp = rtio_input_data(self.ttl_chan) # eat counts of RECENTLY CLOSED sub-bin
            total_counts += counts_tmp
            bin_counter += 1

            '''
            DYNAMIC PROCESSING
            '''
            # handle potential cases where we have yuuug counts
            if counts_tmp <= self.max_counts_bin:
                counts_tmp = self.max_counts_bin
            # update bright/dark likelihoods recursively
            # note: we ignore dark => bright decays for extreme simplicity
            p_b *= self.count_rate_bright_bin[counts_tmp]
            p_d *= self.count_rate_bright_bin[counts_tmp]

            # completion condition - bright state
            if p_d < (p_b * self.error_threshold_frac):
                ion_state = 1
                break
            # completion condition - dark state
            elif p_b < (p_d * self.error_threshold_frac):
                ion_state = 0
                break

        # ensure all gates are closed
        rtio_output(self.ttl_chan, CONFIG_SEND_COUNT_EVENT)
        self.core.break_realtime()

        # eat all remaining counts
        # todo: is this necessary? pretty fair to say it's certain we only have ONE final bin
        count_events_remaining = self.fetch_timestamped_count(now_mu() + 10000)
        while count_events_remaining[0] != -1:
            count_events_remaining = self.fetch_timestamped_count(now_mu() + 10000)
            delay_mu(10000)
        self.core.break_realtime()

        # return results
        return ion_state, total_counts, bin_counter

