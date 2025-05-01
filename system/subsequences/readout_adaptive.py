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

        # timing/binning
        "time_bin_us", "time_bin_mu", "num_sub_bins",

        # count rates & probabilities
        "max_counts_bin", "likelihood_bright_n",
        "likelihood_dark_n", "error_threshold", "error_threshold_frac", "sigma_max", "_t_decay_const"
    }

    def build_subsequence(self, time_bin_us: TFloat = 100, error_threshold: TFloat = 1e-2,
                          sigma_max: TInt32 = 6):
        """
        Defines the main interface for the subsequence.
        Arguments:
            time_bin_us: the individual bin time
            error_threshold: error threshold (fractional) required to determine ion state.
            sigma_max: number of stdevs above mean to account for when processing counts.
        """
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
        Prepare & precompute experimental values.
        """
        '''SANITIZE INPUT'''
        # timing/latencies become challenging below/near 5us bin times
        if self.time_bin_us <= 5: raise ValueError("Invalid bin timing for ReadoutAdaptive. time_bin_us must be >5us.")

        '''PREPARE PARAMETERS'''
        time_readout_us =   self.get_parameter('time_readout_us', group='timing', override=False)
        # rescale count rates for given bin times
        count_rate_bright_bin = (self.get_parameter('count_rate_bright_3ms', group='pmt', override=False) *
                                 (self.time_bin_us * us) / (3. * ms))
        count_rate_dark_bin =   (self.get_parameter('count_rate_dark_3ms', group='pmt', override=False) *
                                 (self.time_bin_us * us) / (3. * ms))

        '''PREPARE HARDWARE'''
        self.ttl_chan_in =  self.pmt.pmt.channel
        self.ttl_chan_out = self.pmt.pmt.channel << 8

        self.time_bin_mu =  self.core.seconds_to_mu(self.time_bin_us * us)
        self.num_sub_bins = time_readout_us // self.time_bin_us

        '''PRECALCULATE STATISTICS'''
        # calculate values to account for Bright => Dark decay
        self._t_decay_const = self.time_bin_us * us / 1.149   # D-5/2 to S-1/2 decay time

        # reprocess error threshold for faster calculation (avoids expensive divisions during kernel)
        self.error_threshold_frac = self.error_threshold / (1. - self.error_threshold)
        # set max counts as given \sigma above distribution parameter
        self.max_counts_bin = round(
            (count_rate_bright_bin + count_rate_dark_bin) +
            self.sigma_max * math.sqrt(count_rate_bright_bin + count_rate_dark_bin)
        )

        # calculate likelihood distributions and store
        likelihood_bright_n, likelihood_dark_n = self._prepare_likelihoods(
            count_rate_bright_bin, count_rate_dark_bin, self.max_counts_bin
        )
        self.set_dataset("likelihood_bright_n", likelihood_bright_n)
        self.setattr_dataset("likelihood_bright_n")
        self.set_dataset("likelihood_dark_n", likelihood_dark_n)
        self.setattr_dataset("likelihood_dark_n")

    @rpc
    def _prepare_likelihoods(self, count_rate_bright: TFloat, count_rate_dark: TFloat,
                             max_counts_bin: TInt32) ->TTuple([TArray(TFloat, 1), TArray(TFloat, 1)]):
        """
        Precalculate the target likelihood distributions for all possible count values
        to avoid expensive on-kernel calculation.
        Arguments:
            count_rate_bright: mean counts per bin for a bright state.
            count_rate_dark: mean counts per bin for a dark state.
            max_counts_bin: maximum number of counts to consider for likelihood calculation.
        Returns:
            array of likelihoods for all possible detected counts (until max_counts_bin) for
            bright and dark states, separately.
        """
        # precalculate factorials & likelihood functions
        # note: use stirling's approx to 2nd order to reduce size of numbers b/c we reach 64b quickly
        log_factorial_stirling = lambda n: (
                math.log(n) * (n + 0.5) - n + 0.5*math.log(2.*math.pi) +
                math.log(1. + 1./(12.*n))
        )
        likelihood_poisson =        lambda ld, n: ld ** n * math.exp(-ld) / math.factorial(n)
        log_likelihood_poisson =    lambda ld, n: n * math.log(ld) - ld - log_factorial_stirling(n)

        # use log-likelihoods to reduce numerical errors b/c numbers are very large BEFORE division
        likelihood_bright_n = np.array([
            likelihood_poisson(count_rate_bright, n) if n <= 10
            else math.exp(log_likelihood_poisson(count_rate_bright, n))
            for n in range(0, max_counts_bin + 1)
        ])
        likelihood_dark_n = np.array([
            likelihood_poisson(count_rate_dark, n) if n <= 10
            else math.exp(log_likelihood_poisson(count_rate_dark, n))
            for n in range(0, max_counts_bin + 1)
        ])
        return likelihood_bright_n, likelihood_dark_n


    '''
    HARDWARE METHODS
    '''
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
        '''PREPARE'''
        # instantiate variables
        ion_state =     -1
        total_counts =  0
        bin_counter =   0
        p_b =           1.  # likelihood bright
        p_d =           1.  # likelihood dark
        # m_n =           1.  # running m (dark)
        # s_n =           0.  # running s (bright => dark)

        # set readout profile and turn on readout beams
        self.pump.readout()
        self.pump.on()
        self.repump_cooling.on()

        # get fiducial start time
        time_start_mu = now_mu()

        # start initial sub-bin (to give slack later on)
        at_mu(time_start_mu)
        rtio_output(self.ttl_chan_out, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)
        delay_mu(self.time_bin_mu)
        rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)

        '''DYNAMIC PROCESSING'''
        while bin_counter < self.num_sub_bins:
            bin_counter += 1

            # schedule next bin (allows extra slack without real overhead)
            at_mu(time_start_mu + bin_counter * (self.time_bin_mu + 8))
            rtio_output(self.ttl_chan_out, CONFIG_COUNT_RISING | CONFIG_RESET_TO_ZERO)
            delay_mu(self.time_bin_mu)
            rtio_output(self.ttl_chan_out, CONFIG_SEND_COUNT_EVENT)

            # eat counts from previous bin and update loop
            counts_tmp = rtio_input_data(self.ttl_chan_in)
            total_counts += counts_tmp
            # handle potential cases where we have yuuug counts
            if counts_tmp >= self.max_counts_bin:
                counts_tmp = self.max_counts_bin

            # update bright/dark likelihoods recursively
            # note: we ignore dark => bright decays for extreme simplicity
            p_b *= self.likelihood_bright_n[counts_tmp] # this is B(n)
            p_d *= self.likelihood_dark_n[counts_tmp]   # this is D(n)

            # # tmp remove - new math
            # # update bright probability
            # p_b *= self.likelihood_bright_n[counts_tmp]
            # # update recursive variables for dark probability
            # s_n = (s_n + m_n) * self.likelihood_bright_n[counts_tmp]
            # m_n *= self.likelihood_dark_n[counts_tmp]
            # p_d = (1. - bin_counter * self._t_decay_const) * m_n + self._t_decay_const * s_n
            # # tmp remove - new math

            # completion condition - bright state
            if p_d < (p_b * self.error_threshold_frac):
                ion_state = 1
                break
            # completion condition - dark state
            elif p_b < (p_d * self.error_threshold_frac):
                ion_state = 0
                break

        # ensure remaining count bin is cleared from input FIFO
        rtio_input_data(self.ttl_chan_in)
        delay_mu(5000)
        return ion_state, total_counts, bin_counter, p_b, p_d

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Ensure all input events are cleared/eaten (in case reset isn't called).
        """
        self.core.break_realtime()

        # eat any remaining count events
        count_events_remaining = self.pmt.fetch_timestamped_count(now_mu() + 10000)
        while count_events_remaining[0] != -1:
            count_events_remaining = self.pmt.fetch_timestamped_count(now_mu() + 10000)
            delay_mu(10000)
        self.core.break_realtime()

