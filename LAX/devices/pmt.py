from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice


class PMTCounter(LAXDevice):
    """
    A wrapper around the TTL edge_counter object for the PMT.
    """

    device_names = {'pmt': 'ttl_counter0'}
    device_parameters = {
        'gating_edge': ('pmt.pmt_gating_edge', None)
    }


    @kernel(flags='fast-math')
    def prepare_devices(self):
        # get default gating edge for counting
        self.counting_method = getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')


    @kernel
    def count(self, time_mu):
        """
        Counts the specified gating edges for a given time.
        """
        self.counting_method(time_mu)
