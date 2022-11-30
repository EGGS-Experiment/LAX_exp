from artiq.experiment import *
from LAX_exp.LAX.base_classes import *


class PMT_counter(object):
    """
    A wrapper around the TTL edge_counter object for the PMT.
    """

    TTL_CHANNEL =       'pmt.pmt_input_channel'
    TTL_GATING_EDGE =   'pmt.pmt_gating_edge'


    # SETUP
    def build(self):
        # get core device
        self.setattr_device("core")

        # get TTL counter channel
        ttl_channel = self.get_dataset(self.TTL_CHANNEL, archive=False)
        self.dev = self.get_device('ttl_counter{:d}'.format(ttl_channel))

        # get default gating edge for counting
        gating_method = self.get_dataset(self.TTL_GATING_EDGE, archive=False)
        self.counting_method = getattr(self.dev, "gate_{:s}_mu".format(gating_method))

    def __getattr__(self, attribute):
        """
        Call methods of the backing TTL counter if not otherwise implemented.
        """
        return getattr(self.dev, attribute)


    # CONVENIENCE METHODS
    @kernel
    def count(self, time_mu):
        """
        # todo: document
        """
        self.counting_method(time_mu)
