from numpy import int32
from artiq.experiment import *


class MeanFilter(HasEnvironment):
    """
    Helper: Mean Filter

    Implements a mean filter (i.e. fixed-length moving average) for user convenience.
    This is separate from FIRFilter for optimization purposes.
    Written to allow operation in an artiq kernel.
    """
    name = 'Mean Filter'
    kernel_invariants = {
        "max_filter_length", "filter_length"
    }

    def build(self, filter_length: TInt32=1) -> TNone:
        """
        Build the Mean filter object.
        :param filter_length: length of the moving average.
        """
        # store build arguments
        self.filter_length = filter_length
        # hardcoded variables
        self.max_filter_length = 1000

    def prepare(self):
        """
        Prepare hardcoded values ahead of time.
        """
        self._prepare_argument_checks()

        # create dynamic filter variables
        self._filter_hist = [int32(0)] * self.filter_length
        self._filter_idx =  int32(0)
        self._filter_out =  int32(0) # note: use int32 for speed since we probably don't need to process large numbers

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure filter_taps meets criteria
        if not (1 <= self.filter_length <= self.max_filter_length):
            raise ValueError("filter_length outside valid range."
                             "Must be in range [1, {:d}].".format(self.max_filter_length))


    '''
    FILTER INTERFACE
    '''
    @portable(flags={"fast-math"})
    def update_single(self, value: TInt32) -> TNone:
        """
        Update the Mean filter with a single value.
        Filter output is calculated immediately upon update.
        :param value: Value to update the filter with.
        """
        # update filter history and running average
        self._filter_hist[self._filter_idx % self.filter_length] = value
        self._filter_out += value
        self._filter_idx += 1

        # subtract history from running average (i.e. circular buffer) once we have stored enough counts
        # note: don't need to check for current index to exceed length since higher values should be all 0
        self._filter_out -= self._filter_hist[self._filter_idx % self.filter_length]

    @portable(flags={"fast-math"})
    def get_current(self) -> TInt32:
        """
        Simply return the current value of the MeanFilter.
        Does not recalculate the filter output.
        :return: The current value of the MeanFilter.
        """
        return self._filter_out

    @portable(flags={"fast-math"})
    def clear(self) -> TNone:
        """
        Reset the Mean filter.
        Clears the running count as well as the storage array.
        """
        self._filter_idx = 0
        self._filter_out = 0
        for i in range(self.filter_length):
            self._filter_hist[i] = 0


if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from time import perf_counter

    '''PREPARE'''
    # create test filter
    _filter_length = 200
    filter_test = MeanFilter((None, None, None, None), filter_length=_filter_length)
    filter_test.prepare()

    # create test sequence
    test_input = np.concatenate([np.ones(1000),
                                 np.zeros(250),
                                 np.ones(50)])
    test_output = np.zeros(len(test_input))

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

        # process output - todo
        # if (i > _filter_length) and (status_now != status_latched):
        #     filter_test.clear()


    '''ANALYZE'''
    # print timing overheads
    print("\n\tTotal Overhead:\t\t{:.3f} us".format(overhead_total_mu * 1e6))
    print("\tAverage Overhead:\t{:.3f} ns\n".format(overhead_total_mu / len(test_input) * 1e9))

    # plot results
    plt.plot(test_input, label='Input')
    plt.plot(test_output / _filter_length, label='Output')
    plt.legend(loc='best')
    plt.show()

