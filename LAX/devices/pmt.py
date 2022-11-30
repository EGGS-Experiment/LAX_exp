from artiq.experiment import *
from inspect import getmembers, ismethod


class PMT_counter(HasEnvironment):
    """
    A wrapper around the TTL edge_counter object for the PMT.
    """

    kernel_invariants =     set()

    TTL_CHANNEL =       'pmt.pmt_input_channel'
    TTL_GATING_EDGE =   'pmt.pmt_gating_edge'


    # SETUP
    def build(self):
        # get core device
        self.setattr_device("core")

        # get TTL counter channel
        # ttl_channel = self.get_dataset(self.TTL_CHANNEL, archive=False)
        # setattr(
        #     self,
        #     "dev",
        #     self.get_device('ttl_counter{:d}'.format(ttl_channel))
        # )
        #
        # # get default gating edge for counting
        # gating_method = self.get_dataset(self.TTL_GATING_EDGE, archive=False)
        # setattr(
        #     self,
        #     "counting_method",
        #     getattr(self.dev, "gate_{:s}_mu".format(gating_method))
        # )
        ttl_channel = self.get_dataset(self.TTL_CHANNEL, archive=False)
        self.dev = self.get_device('ttl_counter{:d}'.format(ttl_channel))

        # get default gating edge for counting
        gating_method = self.get_dataset(self.TTL_GATING_EDGE, archive=False)
        self.counting_method = getattr(self.dev, "gate_{:s}_mu".format(gating_method))

        # steal all relevant methods of underlying TTL counter object so users
        # can directly call methods from this wrapper
        isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")
        device_functions = getmembers(self.dev, isDeviceFunction)
        for (function_name, function_object) in device_functions:
            setattr(self, function_name, function_object)


    # CONVENIENCE METHODS
    @kernel
    def count(self, time_mu):
        """
        # todo: document
        """
        self.counting_method(time_mu)
