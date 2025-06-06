from artiq.experiment import *
from numpy import ndarray, int32, int64
# todo: handle floats somehow, or stick with int32? idk


class FIRFilter(HasEnvironment):
    """
    Helper: FIR Filter

    Implements an FIR filter for user convenience.
    Written to allow operation in an artiq kernel.
    """
    name = 'FIR Filter'
    kernel_invariants = {
        "filter_taps", "max_filter_length", "_filter_length"
    }

    def build(self, filter_taps) -> TNone:
        """
        Build the FIR filter object.
        :param filter_taps: a list of the filter taps for the FIR filter.
        """
        # get relevant devices
        # todo: is this necessary?
        # self.setattr_device('core')

        # store build arguments
        self.filter_taps = filter_taps

        # hardcoded variables
        self.max_filter_length = 1000

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # # ensure filter_taps is of correct type
        # if not isinstance(self.filter_taps, list) or not isinstance(self.filter_taps, ndarray):
        #     raise ValueError("Invalid type for filter_taps. filter_taps must be a list object.")
        # if isinstance(self.filter_taps, list) and not all(isinstance(x, int32) for x in self.filter_taps):
        #     raise ValueError("Invalid data types in filter_taps. All values in filter_taps must be of type int32.")
        # if isinstance(self.filter_taps, ndarray) and not (self.filter_taps.dtype == int32):
        #     raise ValueError("Invalid data types in filter_taps. All values in filter_taps must be of type int32.")

        # ensure filter_taps meets criteria
        if not (1 <= len(self.filter_taps) <= self.max_filter_length):
            raise ValueError("filter_taps outside valid range. Must be in range [1, {:d}].".format(self.max_filter_length))

    def prepare(self):
        """
        Prepare hardcoded values ahead of time.
        """
        self._prepare_argument_checks()

        # hardcoded variables
        self._filter_length = len(self.filter_taps)
        # reverse filter taps to speed up calculation
        self.filter_taps.reverse()

        # dynamic filter variables
        self._filter_hist = [int32(0)] * self._filter_length # stores history of input values
        self._filter_idx =  int32(0)    # stores running count of the filter
        self._filter_updated = False    # flag to indicate whether output needs to be recalculated
        # note: use int32 for _filter_out for speed since we probably don't need to process large numbers
        self._filter_out = int32(0)     # stores the current output of the FIRFilter


    '''
    FILTER INTERFACE
    '''
    @portable(flags={"fast-math"})
    def update_single(self, value: TInt32) -> TNone:
        """
        Update the FIR filter with a single value.
        To minimize overhead, this function does not calculate the filter output.
        :param value: Value to update the filter with.
        """
        self._filter_hist[self._filter_idx % self._filter_length] = value
        self._filter_updated = True
        self._filter_idx += 1

    @portable(flags={"fast-math"})
    def get_current(self) -> TInt32:
        """
        Return current output of the FIRFilter.
        To minimize overheads, filter output is only calculated on the first call
            to this function (until the filter is updated); successive calls to this function
            will simply return the previously calculated value.
        :return: The current FIRFilter output value.
        """
        # recalculate filter output only if we have new data
        if self._filter_updated is True:
            # calculate convolution
            self._filter_out = 0
            for i in range(self._filter_length):
                self._filter_out += (self._filter_hist[i] *
                                     self.filter_taps[(i + self._filter_idx - 1) % self._filter_length])

            # set update flag
            self._filter_updated = False

        return self._filter_out

    @portable(flags={"fast-math"})
    def clear(self) -> TNone:
        """
        Reset the FIR filter.
        Clears the running count as well as the storage array.
        """
        self._filter_idx = 0
        self._filter_out = 0
        self._filter_updated = False

        for i in range(self._filter_length):
            self._filter_hist[i] = 0


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from time import perf_counter

    '''PREPARE'''
    # user configuration
    f_sample = 44100
    time_length = 0.01
    freq_arrs = [1e2, 1e3, 1e4, 1.5e4, 2e4]

    # create FIR filter sequence
    filter_taps = [
        0.004087, 0.001930,-0.001847,-0.005502,-0.006956,0.004898,0.000307,0.006469,0.010407,0.009493,
        0.003213,-0.006136,-0.014116,-0.016107,-0.009740,0.003374,0.017666,0.025672,0.021492,0.004177,
        -0.020612,-0.041731,-0.046388,-0.025354,0.022561,0.088761,0.156853,0.207858,0.226757,0.207858,
        0.156853,0.088761,0.022561,-0.025354,-0.046388,-0.041731,-0.020612,0.004177,0.021492,0.025672,
        0.017666,0.003374,-0.009740,-0.016107,-0.014116,-0.006136,0.003213,0.009493,0.010407,0.006469,
        0.000307,-0.004898,-0.006956,-0.005502,-0.001847,0.001930,0.004087
    ]   # lpf w/cutoff @ 0.11337 * f_sample (5 kHz for f_sample @ 44.1 kHz)

    # create test filter object
    filter_test = FIRFilter((None, None, None, None), filter_taps=filter_taps)
    filter_test.prepare()

    # create test sequence
    t_arr = np.arange(0., time_length, 1. / f_sample)
    num_samples = len(t_arr)
    test_input = np.sum([
        np.sin((2.*np.pi*f_hz) * (t_arr + np.random.normal(0., 1./f_sample/25., num_samples)))
        for f_hz in freq_arrs
    ], axis=0)
    test_output = np.zeros(num_samples)

    # prepare to benchmark timing overheads/performance
    time_start_mu, overhead_total_mu = (0., 0.)


    '''RUN'''
    # apply filter on test sequence
    (status_now, status_latched) = (0, 0)
    for i, val in enumerate(test_input):
        time_start_mu = perf_counter()
        filter_test.update_single(val)
        filter_out = filter_test.get_current()
        overhead_total_mu += perf_counter() - time_start_mu

        test_output[i] = filter_out


    '''ANALYZE'''
    # print timing overheads
    print("\n\tTotal Overhead:\t\t{:.3f}\tms".format(overhead_total_mu * 1e3))
    print("\tAverage Overhead:\t{:.3f}\tus\n".format(overhead_total_mu / num_samples * 1e6))

    # plot time-domain results
    plt.plot(test_input, label='Input')
    plt.plot(test_output, label='Output')
    plt.title("Time Domain")
    plt.legend(loc='best')
    plt.show()

    # plot frequency domain results
    fft_input = np.abs(np.fft.fft(test_input)[:num_samples//2])
    fft_output = np.abs(np.fft.fft(test_output)[:num_samples//2])
    plt.plot(fft_input, label='Input')
    plt.plot(fft_output, label='Output')
    plt.title("Frequency Domain")
    plt.yscale('log')
    plt.legend(loc='best')
    plt.show()

