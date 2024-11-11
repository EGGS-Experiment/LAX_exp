import numpy as np
from artiq.experiment import *


class RPCTimeProfiler(EnvExperiment):
    """
    RPC Time Profiling

    Profiles the amount of time taken for an RPC call.
    """

    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("led1")

        # get user arguments
        self.setattr_argument("repetitions",    NumberValue(default=10000, ndecimals=0, step=1000, min=10, max=10000000))
        self.setattr_argument("time_delay_ms",  NumberValue(default=3, ndecimals=6, step=1, min=1, max=100))

    def prepare(self):
        # convert relevant arguments to machine units
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_ms * ms)

        # store machine unit conversion
        self.set_dataset("core_mu_to_seconds", self.core.ref_period)

        # create dataset for storing results
        self.set_dataset("results", np.zeros(self.repetitions, dtype=np.int64))
        self.setattr_dataset("results")
        self._idx_results = 0

    @rpc(flags={"async"})
    def update_results(self, time_mu: TInt64) -> TNone:
        """
        Update the result dataset.
        :param time_mu: the response time, in machine units.
        """
        self.mutate_dataset('results', self._idx_results, time_mu)
        self.set_dataset('management.dynamic.completion_pct',
                         round(self._idx_results / self.repetitions * 100., 3),
                         broadcast=True, persist=True, archive=False)
        self._idx_results += 1

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare system
        self.core.reset()
        self.led1.off()

        # main loop
        for i in range(self.repetitions):
            self.core.break_realtime()

            # simulate pulse sequence
            self.led1.pulse_mu(self.time_delay_mu)

            # RPC call
            time_start_mu = self.core.get_rtio_counter_mu()
            self.rpc_call()
            time_stop_mu = self.core.get_rtio_counter_mu()
            self.update_results(time_stop_mu - time_start_mu)

        # wait until all submitted events complete
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @rpc
    def rpc_call(self) -> TNone:
        """
        Empty RPC call to measure bare overhead of synchronous RPC.
        """
        return

    def analyze(self):
        # convert results from mu to seconds
        results_s = np.array([self.core.mu_to_seconds(time_mu) for time_mu in self.results])

        # process result statistics
        res_mean, res_stdev = (np.mean(results_s), np.std(results_s))
        res_stderr = res_stdev / np.sqrt(len(results_s))
        self.set_dataset("rpc_time_s_mean", res_mean)
        self.set_dataset("rpc_time_s_stdev", res_stdev)
        self.set_dataset("rpc_time_s_stderr", res_stderr)

        # print statistics to log
        print("Results:")
        print("\tRPC Time (ms):\t{:.3f} +/- {:.3f}".format(res_mean / ms, res_stderr / ms))
        print("\tRPC stdev (ms):\t{:.3f}".format(res_stdev / ms))
