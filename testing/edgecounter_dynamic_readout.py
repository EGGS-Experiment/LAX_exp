import math
import numpy as np
from artiq.experiment import *
from artiq.coredevice.rtio import rtio_output, rtio_input_timestamp, rtio_input_data
from artiq.coredevice.edge_counter import (CONFIG_COUNT_RISING, CONFIG_COUNT_FALLING,
                                           CONFIG_SEND_COUNT_EVENT, CONFIG_RESET_TO_ZERO)

# todo: subtract 3x slack time from bin time to make it equal
class TTLDynamicTest(EnvExperiment):
    """
    TTLDynamicTest
    idk 2 TTLDynamicTest
    """
    kernel_invariants = {
        # devices
        "ttl", "ttl_chan_out", "ttl_chan_in",

        # timing
        "time_readout_us", "time_readout_mu", "time_bin_us", "time_bin_mu", "time_bin_process_slack_mu",

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
        self.setattr_device("scheduler")
        self.setattr_device("ttl0_counter")

        self.repetitions =      10000
        self.time_readout_us =  1500
        self.time_bin_us =      100

        self.count_rate_bright =    120     # per 3ms
        self.count_rate_dark =      36      # per 3ms
        self.error_threshold =      1e-2    # total error fraction
        self._dist_sigma_max =      6       # \sigma above mean for max counts

    def prepare(self):
        """
        Prepare & precompute experimental values.
        """
        '''PREPARE HARDWARE'''
        self.ttl =      self.get_device("ttl0_counter")
        self.ttl_chan_in =  self.ttl.channel
        self.ttl_chan_out = self.ttl.channel << 8
        self.time_bin_process_slack_mu = self.core.seconds_to_mu(0.1 * us)

        self.time_bin_mu =      self.core.seconds_to_mu(self.time_bin_us * us)
        self.num_sub_bins =     self.time_readout_us // self.time_bin_us
        self.time_readout_mu =  self.core.seconds_to_mu((self.time_bin_us * us) * self.num_sub_bins)

        '''PRECALCULATE LIKELIHOOD DISTRIBUTIONS'''
        # rescale the count rates for the given bin times
        self.count_rate_bright_bin =    self.count_rate_bright * ((self.time_bin_us * us) / (3. * ms))
        self.count_rate_dark_bin =      self.count_rate_dark * ((self.time_bin_us * us) / (3. * ms))

        # reprocess error threshold for faster calculation (avoids expensive divisions during kernel)
        self.error_threshold_frac =     self.error_threshold / (1. - self.error_threshold)

        self._prepare_likelihoods()
        self._prepare_dataset()

    def _prepare_likelihoods(self):
        """
        Precalculate the target likelihood distributions for all possible count values
        to avoid expensive on-kernel calculation.
        """
        # set max counts as given \sigma above distribution parameter
        self.max_counts_bin = round(
            (self.count_rate_bright_bin + self.count_rate_dark_bin) +
            self._dist_sigma_max * math.sqrt(self.count_rate_bright_bin + self.count_rate_dark_bin)
        )

        # precalculate factorials & likelihood functions
        log_factorial_stirling = lambda n: (
                math.log(n) * (n + 0.5) - n + 0.5*math.log(2.*math.pi) +
                math.log(1. + 1./(12.*n))
        )
        likelihood_poisson =        lambda ld, n: ld ** n * math.exp(-ld) / math.factorial(n)
        log_likelihood_poisson =    lambda ld, n: n * math.log(ld) - ld - log_factorial_stirling(n)

        # use log-likelihoods b/c numbers are very large BEFORE division
        # this reduces numerical digitization errors
        self.likelihood_bright_n = np.array([
            likelihood_poisson(self.count_rate_bright_bin, n) if n <= 10
            else math.exp(log_likelihood_poisson(self.count_rate_bright_bin, n))
            for n in range(0, self.max_counts_bin + 1)
        ])
        self.likelihood_dark_n = np.array([
            likelihood_poisson(self.count_rate_dark_bin, n) if n <= 10
            else math.exp(log_likelihood_poisson(self.count_rate_dark_bin, n))
            for n in range(0, self.max_counts_bin + 1)
        ])
        self.set_dataset("likelihood_bright_n", self.likelihood_bright_n)
        self.setattr_dataset("likelihood_bright_n")
        self.set_dataset("likelihood_dark_n", self.likelihood_dark_n)
        self.setattr_dataset("likelihood_dark_n")

    def _prepare_dataset(self):
        """
        Prepare datasets to record data.
        """
        # create data structures for results
        self.set_dataset('results', np.zeros((self.repetitions, 5)))
        self.setattr_dataset('results')
        self._result_iter = 0
        self._counts_iter = 0

        # create data structures for dynamic updates
        self._dynamic_reduction_factor = self.get_dataset('management.dynamic_plot_reduction_factor',
                                                          default=10, archive=False)
        self.kernel_invariants.add("_dynamic_reduction_factor")
        self._completion_iter_to_pct = 100. / len(self.results)
        self.kernel_invariants.add("_completion_iter_to_pct")

        # set up dynamic datasets
        self.set_dataset('management.dynamic.completion_pct', 0.,
                         broadcast=True, persist=True, archive=False)
        # downsample counts for dynamic plotting
        dynamic_counts_len = (self.repetitions // self._dynamic_reduction_factor) + 1
        dynamic_counts_arr = np.zeros(dynamic_counts_len, dtype=np.int32) * np.nan
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        dynamic_counts_arr[0] = 0
        self.set_dataset('temp.counts.trace', dynamic_counts_arr,
                         broadcast=True, persist=False, archive=False)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset
        self.core.reset()

        # todo: add error handling
        for i in range(self.repetitions):
            # simulate pulse sequence
            self.core.break_realtime()
            delay_mu(1000000)   # 1ms

            # dynamic readout
            results = self._readout()
            self.core.break_realtime()

            # add slack/finish up
            self.update_results(results)
            self.core.break_realtime()

            # periodically check termination
            if i % 10 == 0:
                if self.scheduler.check_termination():
                    self.core.break_realtime()
                    break

        # clean up
        self.core.wait_until_mu(now_mu())
        delay_mu(1000000)
        self.core.break_realtime()
        self.core.reset()


    @rpc(flags={"async"})
    def update_results(self, args):
        # # tmp remove
        # if self._result_iter % 5 == 0:
        #     print(args)
        # # tmp remove

        # store results in main dataset
        self.mutate_dataset('results', self._result_iter, np.array(args))

        # do intermediate processing
        if (self._result_iter % self._dynamic_reduction_factor) == 0:
            # plot counts in real-time to monitor ion death
            self.mutate_dataset('temp.counts.trace', self._counts_iter, args[1])
            self._counts_iter += 1
            # monitor completion status
            self.set_dataset('management.dynamic.completion_pct',
                             round(self._result_iter * self._completion_iter_to_pct, 3),
                             broadcast=True, persist=True, archive=False)

        # increment result iterator
        self._result_iter += 1

    @kernel(flags={"fast-math"})
    def _readout(self) -> TTuple([TInt32, TInt32, TInt32, TFloat, TFloat]):
        """
        Dynamically reads out the ion state using an adaptive maximum likelihood technique.
        See A.H. Burrell's thesis (Oxford, 2010) for more details.
        Returns: a tuple of (ion_state, total_counts, elapsed_bins, P_bright, P_dark).
        ion_state: the determined ion state - can be any of (-1, 0, 1). -1 means "indeterminate," i.e. the
            bright/dark discrimination criteria was not met. 0 means "dark." 1 means "bright."
        total_counts: the total number of counts detected during readout. This number will depend on the number
            of bins that were required to determine the ion state.
        elapsed_bins: the number of constituent bins required to determine the ion state.
        P_bright: the likelihood that the ion was bright, given the count "trajectory" detected.
        P_bright: the likelihood that the ion was dark, given the count "trajectory" detected.
        """
        # get reference time
        time_start_mu = now_mu()

        # instantiate variables
        ion_state =     -1
        total_counts =  0
        bin_counter =   0
        p_b =           1.  # likelihood bright
        p_d =           1.  # likelihood dark

        # start initial sub-bin
        at_mu(time_start_mu)
        rtio_output(self.ttl_chan_out, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)
        delay_mu(self.time_bin_mu)
        rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)

        # dynamically process each sub-bin
        while bin_counter < self.num_sub_bins:
            # schedule next bin (allows extra slack without real overhead)
            at_mu(time_start_mu + (bin_counter + 1) * (self.time_bin_mu + 8))
            rtio_output(self.ttl_chan_out, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)
            delay_mu(self.time_bin_mu)
            rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)

            # eat counts and update loop
            counts_tmp = rtio_input_data(self.ttl_chan_in) # eat counts of RECENTLY CLOSED sub-bin
            delay_mu(self.time_bin_process_slack_mu)  # tmp remove - add extra slack
            total_counts += counts_tmp
            bin_counter += 1

            '''
            DYNAMIC PROCESSING
            '''
            # handle potential cases where we have yuuug counts
            if counts_tmp >= self.max_counts_bin:
                counts_tmp = self.max_counts_bin
            # update bright/dark likelihoods recursively
            # note: we ignore dark => bright decays for extreme simplicity
            p_b *= self.likelihood_bright_n[counts_tmp]
            p_d *= self.likelihood_dark_n[counts_tmp]

            # completion condition - bright state
            if p_d < (p_b * self.error_threshold_frac):
                ion_state = 1
                break
            # completion condition - dark state
            elif p_b < (p_d * self.error_threshold_frac):
                ion_state = 0
                break

        # ensure all gates are closed
        rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)
        self.core.break_realtime()

        # eat all remaining counts
        # todo: is this necessary? pretty fair to say it's certain we only have ONE final bin
        count_events_remaining = self.ttl.fetch_timestamped_count(now_mu() + 5000)
        while count_events_remaining[0] != -1:
            count_events_remaining = self.ttl.fetch_timestamped_count(now_mu() + 5000)
            delay_mu(5000)
        self.core.break_realtime()

        # return results
        return ion_state, total_counts, bin_counter, p_b, p_d

    def analyze(self):
        """
        Print summary statistics.
        """
        print("\n########## RESULT SUMMARY ##########")

        # collate data
        data_bright =   self.results[self.results[:, 0] == 1]
        data_dark =     self.results[self.results[:, 0] == 0]
        data_idk =      self.results[self.results[:, 0] == -1]

        # print bright/dark/idk percentages
        print("Discrimination Results:\n\tBright:\t{:.6g}%\n\tDark:\t{:.6g}%\n\tidk:\t{:.6g}%".format(
            len(data_bright[:, 0]) / len(self.results) * 100.,
            len(data_dark[:, 0]) / len(self.results) * 100.,
            len(data_idk[:, 0]) / len(self.results) * 100.
        ))

        # print mean count rates
        print("Mean Count Rates (per 3ms):\n\tBright:\t{:.4g}\n\tDark:\t{:.4g}\n\tidk:\t{:.4g}".format(
            np.mean(data_bright[:, 1] / data_bright[:, 2] * (3e-3 / (self.time_bin_us * us))),
            np.mean(data_dark[:, 1] / data_dark[:, 2] * (3e-3 / (self.time_bin_us * us))),
            np.mean(data_idk[:, 1] / data_idk[:, 2] * (3e-3 / (self.time_bin_us * us)))
        ))

        # print mean count rates
        print("Mean Time to Detection (# bins, us):"
              "\n\tBright:\t{:.4g}/\t{:.4g}\n\tDark:\t{:.4g}/\t{:.4g}\n\tidk:\t{:.4g}/\t{:.4g}".format(
            np.mean(data_bright[:, 2]), np.mean(data_bright[:, 2]) * self.time_bin_us,
            np.mean(data_dark[:, 2]), np.mean(data_dark[:, 2]) * self.time_bin_us,
            np.mean(data_idk[:, 2]), np.mean(data_idk[:, 2]) * self.time_bin_us
        ))

