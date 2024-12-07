from artiq.experiment import *
from LAX_exp.base import LAXSubsequence


class NoOperation(LAXSubsequence):
    """
    Subsequence: No Operation

    A "no_op" subsequence to use as a placeholder or null (i.e. do-nothing) sequence.
    """
    name = 'no_operation'

    def build_subsequence(self):
        pass

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        pass
