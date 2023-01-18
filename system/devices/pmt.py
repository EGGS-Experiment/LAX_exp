from artiq.experiment import *
from LAX_exp.base import LAXDevice


class PMTCounter(LAXDevice):
    """
    Device: PMT photon counter

    Wrapper for the TTL edge_counter object for the PMT.
    """
    name = "pmt"

    parameters = {
        'gating_edge':          ('pmt.gating_edge', None)
    }
    core_devices = {
        'pmt': 'ttl0_counter'
    }

    def prepare_device(self):
        self.gating_edge = self.get_parameter('gating_edge', group='pmt', override=False)

        # get default gating edge for counting
        self.counting_method = getattr(self.pmt, 'gate_{:s}_mu'.format(self.gating_edge))
        self.kernel_invariants.add('counting_method')

    @kernel
    def count(self, time_mu):
        """
        Counts the specified gating edges for a given time.
        """
        self.counting_method(time_mu)
