from artiq.experiment import *
from LAX_exp.base import LAXDevice


class ADC(LAXDevice):
    """
    Device: ADC (Analog-Digital Converter)

    Wrapper for the Sampler object.
    """
    name = "adc"
    core_device = ('adc', 'sampler0')

    def prepare_device(self):
        self.gating_edge =                              self.get_parameter('gating_edge', group='pmt', override=False)

        # get default gating edge for counting
        self.counting_method =                          getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')

    @kernel(flags={"fast-math"})
    def count(self, time_mu):
        """
        Counts the specified gating edges for a given time.
        """
        self.counting_method(time_mu)
