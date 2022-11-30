from artiq.experiment import *
import inspect

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

        # get methods
        ppsh41 = lambda okth: (callable(okth)) and (inspect.ismethod(okth)) and (okth.__name__ is not "__init__")
        kka = inspect.getmembers(self.dev, ppsh41)
        for (d_name, d_meth) in kka:
            print('\t{}: {}'.format(d_name, d_meth.__func__))
            setattr(self, d_name, d_meth)

    def __getattr__(self, attr):
        """
        Call methods of the backing TTL counter if not otherwise implemented.
        """
        @kernel
        def methodtmp(th1):
            yz1 = getattr(self.dev, attr)
            return yz1(th1)

        return methodtmp



    # CONVENIENCE METHODS
    @kernel
    def count(self, time_mu):
        """
        # todo: document
        """
        self.counting_method(time_mu)
