from artiq.experiment import *
from LAX_exp.base import LAXExperiment
# change name of subsequence to avoid clashing with our cleanup utility class
from LAX_exp.system.subsequences import Cleanup as CleanupSubsequence


class Cleanup(LAXExperiment, Experiment):
    """
    Utility: Cleanup

    Clean up and reset all hardware in the event of some experiment failure.
    """
    name = 'Cleanup'


    def build_experiment(self):
        self.cleanup_subsequence = CleanupSubsequence(self)

    @property
    def results_shape(self):
        return (1, 1)

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()
        self.core.break_realtime()
        self.cleanup_subsequence.run()
