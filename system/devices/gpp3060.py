from LAX_exp.base import LAXDevice
from artiq.experiment import *

from os import environ
import labrad


class GPP3060(LAXDevice):
    """
    High-level api functions for using the Oven
    """

    name = "gpp3060"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.gpp3060 = self.cxn.GPP3060_server

        self.ovenChannel = 1

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @rpc
    def turn_oven_on(self) -> TNone:
        """
        Turn on the oven to 1 volt
        """
        self.gpp3060.channelVoltage(self.ovenChannel, 1.)

    @rpc
    def turn_oven_off(self) -> None:
        """
        Turn off the oven
        """
        self.gpp3060.channelVoltage(self.ovenChannel, 0)

    @rpc
    def set_oven_voltage(self, voltage):

        """
        Set the oven voltage
        """
        # assert 0. <= voltage <= 30., f"voltage must be between 0V and 30V"
        self.gpp3060.channelVoltage(self.ovenChannel, voltage)

    @rpc
    def set_oven_current(self, current):

        """
        Set the oven voltage
        """
        # assert 0. <= current <= 6., f"current must be between 0A and 6A"
        self.gpp3060.channelCurrent(self.ovenChannel, current)
