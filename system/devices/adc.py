from artiq.experiment import *
from LAX_exp.base import LAXDevice


class ADC(LAXDevice):
    """
    Device: ADC (Analog-Digital Converter)

    Wrapper for the Sampler object.
    """
    name = "adc"
    core_device = ('adc', 'sampler0')
    kernel_invariants = {
        "gating_edge",
        "counting_method"
    }

    def prepare_device(self):
        self.gating_edge = self.get_parameter('gating_edge', group='devices.pmt', override=False)

        # get default gating edge for counting
        self.counting_method = getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))

    @kernel(flags={"fast-math"})
    def count(self, time_mu: TInt64) -> TNone:
        """
        Counts the specified gating edges for a given time.
        Arguments:
            time_mu (TInt64): the time to count for.
        """
        self.counting_method(time_mu)
