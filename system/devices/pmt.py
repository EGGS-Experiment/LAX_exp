from artiq.experiment import *
from LAX_exp.base import LAXDevice


class PMTCounter(LAXDevice):
    """
    A wrapper around the TTL edge_counter object for the PMT.
    """
    name = "pmt"

    parameters = {
        'gating_edge':          ('pmt.pmt_gating_edge', None)
    }
    core_devices = {
        'pmt': 'ttl_counter0'
    }

    def prepare_device(self):
        # get default gating edge for counting
        self.counting_method = getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')

    @kernel
    def count(self, time_mu):
        """
        Counts the specified gating edges for a given time.
        """
        self.counting_method(time_mu)
