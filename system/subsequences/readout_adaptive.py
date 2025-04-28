import math
import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from artiq.coredevice.rtio import rtio_output, rtio_input_data
from artiq.coredevice.edge_counter import (CONFIG_COUNT_RISING, CONFIG_COUNT_FALLING,
                                           CONFIG_SEND_COUNT_EVENT, CONFIG_RESET_TO_ZERO)


class ReadoutAdaptive(LAXSubsequence):
    """
    Subsequence: Readout (Adaptive)

    Dynamically reads out the ion state via state-selective, laser-induced fluorescence.
    Implements an adaptive maximum likelihood technique to reduce readout times.
    See A.H. Burrell's thesis (Oxford, 2010) for more details.
    """
    name = 'readout'
    kernel_invariants = {
        # devices
        "ttl_chan_out", "ttl_chan_in",

        # timing
        "time_readout_mu", "time_bin_us", "time_bin_mu", "time_bin_process_slack_mu",

        # config/binning
        "num_sub_bins",

        # count rates & probabilities
        "count_rate_bright_bin", "count_rate_dark_bin", "max_counts_bin", "likelihood_bright_n",
        "likelihood_dark_n", "error_threshold", "error_threshold_frac", "sigma_max"
    }

    def build_subsequence(self, time_bin_us: TFloat = 100, error_threshold: TFloat = 1e-2,
                          sigma_max: TInt32 = 6):
        """
        Defines the main interface for the subsequence.
        Arguments:
            time_bin_us: the individual bin time
            error_threshold: error threshold (fractional) to determine ion state.
            sigma_max: number of stdevs above mean to account for when processing counts.
        """
        # check argument validity
        if time_bin_us <= 5:
            raise ValueError("Invalid bin timing for ReadoutAdaptive. time_bin_us must be >5us.")

        # set subsequence arguments
        self.time_bin_us =      time_bin_us
        self.error_threshold =  error_threshold
        self.sigma_max =        sigma_max

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

    def prepare_subsequence(self):
        """
        todo document
        """
        '''PREPARE PARAMETERS'''
        # get parameters from dataset manager
        time_readout_us =   self.get_parameter('time_readout_us', group='timing', override=False)
        count_rate_bright = self.get_parameter('count_rate_bright_3ms', group='pmt', override=False)
        count_rate_dark =   self.get_parameter('count_rate_dark_3ms', group='pmt', override=False)

        '''PREPARE HARDWARE'''
        # get RTIO addresses for PMT so we can access it directly
        self.ttl_chan_in =  self.pmt.channel
        self.ttl_chan_out = self.pmt.channel << 8
        self.time_bin_process_slack_mu = self.core.seconds_to_mu(0.5 * us)

        self.time_bin_mu =      self.core.seconds_to_mu(self.time_bin_us * us)
        self.num_sub_bins =     time_readout_us // self.time_bin_us
        self.time_readout_mu =  self.core.seconds_to_mu((self.time_bin_us * us) * self.num_sub_bins)

        '''PRECALCULATE LIKELIHOOD DISTRIBUTIONS'''
        # subtract dark counts from bright to get actual ION SIGNAL, and rescale count rates for specified bin times
        # (we let users specify total bright counts (i.e. ion signal + background) for convenience)
        self.count_rate_bright_bin =    (count_rate_bright - count_rate_dark) * ((self.time_bin_us * us) / (3. * ms))
        self.count_rate_dark_bin =      count_rate_dark * ((self.time_bin_us * us) / (3. * ms))

        # reprocess error threshold for faster calculation (avoids expensive divisions during kernel)
        self.error_threshold_frac =     self.error_threshold / (1. - self.error_threshold)

        # calculate likelihood statistics for each possible bin count
        self._prepare_likelihoods()

    def _prepare_likelihoods(self):
        """
        Precalculate the target likelihood distributions for all possible count values
        to avoid expensive on-kernel calculation.
        """
        # set max counts as given \sigma above distribution parameter
        self.max_counts_bin = round(
            (self.count_rate_bright_bin + self.count_rate_dark_bin) +
            self.sigma_max * math.sqrt(self.count_rate_bright_bin + self.count_rate_dark_bin)
        )

        # precalculate factorials & likelihood functions
        log_factorial_stirling = lambda n: (
                math.log(n) * (n + 0.5) - n + 0.5*math.log(2.*math.pi) +
                math.log(1. + 1./(12.*n))
        )
        likelihood_poisson =        lambda ld, n: ld ** n * math.exp(-ld) / math.factorial(n)
        log_likelihood_poisson =    lambda ld, n: n * math.log(ld) - ld - log_factorial_stirling(n)

        # use log-likelihoods to reduce (numerical digitization errors) b/c numbers are VERY large BEFORE division
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

    @kernel(flags={"fast-math"})
    def run(self) -> TTuple([TInt32, TInt32, TInt32, TFloat, TFloat]):
        """
        Dynamically reads out the ion state using an adaptive maximum likelihood technique.
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
        ion_state =     -1  # ion state starts unknown
        total_counts =  0   # total number of collected counts
        bin_counter =   0   # elapsed bins
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

            # check completion condition - bright state
            if p_d < (p_b * self.error_threshold_frac):
                ion_state = 1
                break
            # check completion condition - dark state
            elif p_b < (p_d * self.error_threshold_frac):
                ion_state = 0
                break

        # ensure all gates are closed and eat remaining count bin
        rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)
        rtio_input_data(self.ttl_chan_in)
        # todo: check if this break_realtime() is necessary, or if it can be replaced w/fixed slack
        # self.core.break_realtime()

        # return results
        return ion_state, total_counts, bin_counter, p_b, p_d

    def analyze_subsequence(self):
        """
        Print summary statistics.
        # todo - need to store results on our own if we want to do statistics on them
        """
        pass
        # print("\n########## RESULT SUMMARY ##########")
        #
        # # collate data
        # data_bright =   self.results[self.results[:, 0] == 1]
        # data_dark =     self.results[self.results[:, 0] == 0]
        # data_idk =      self.results[self.results[:, 0] == -1]
        #
        # # todo: add stds for calculating these
        #
        # # print bright/dark/idk percentages
        # print("Discrimination Results:\n\tBright:\t{:.6g}%\n\tDark:\t{:.6g}%\n\tidk:\t{:.6g}%".format(
        #     len(data_bright[:, 0]) / len(self.results) * 100.,
        #     len(data_dark[:, 0]) / len(self.results) * 100.,
        #     len(data_idk[:, 0]) / len(self.results) * 100.
        # ))
        #
        # # print mean count rates
        # print("Count Rates (per 3ms):\n\tBright:\t{:.4g}\n\tDark:\t{:.4g}\n\tidk:\t{:.4g}".format(
        #     np.mean(data_bright[:, 1] / data_bright[:, 2] * (3e-3 / (self.time_bin_us * us))),
        #     np.mean(data_dark[:, 1] / data_dark[:, 2] * (3e-3 / (self.time_bin_us * us))),
        #     np.mean(data_idk[:, 1] / data_idk[:, 2] * (3e-3 / (self.time_bin_us * us)))
        # ))
        #
        # # print mean count rates
        # print("Time to Detection (# bins, us):"
        #       "\n\tBright:\t{:.4g}/\t{:.4g}\n\tDark:\t{:.4g}/\t{:.4g}\n\tidk:\t{:.4g}/\t{:.4g}".format(
        #     np.mean(data_bright[:, 2]), np.mean(data_bright[:, 2]) * self.time_bin_us,
        #     np.mean(data_dark[:, 2]), np.mean(data_dark[:, 2]) * self.time_bin_us,
        #     np.mean(data_idk[:, 2]), np.mean(data_idk[:, 2]) * self.time_bin_us
        # ))
