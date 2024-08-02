from LAX_exp.base import LAXDevice
from artiq.experiment import *

from os import environ
import labrad


class Oven(LAXDevice):
    """
    High-level api functions for using the Oven
    """

    name = "oven"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.oven = self.cxn.gpp3060_server

        self.ovenChannel = 1

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        pass

    @rpc
    def on(self) -> TNone:
        """
        Toggle Oven off
        """
        self.oven.channel_toggle(self.ovenChannel, 1)

    @rpc
    def off(self) -> None:
        """
        Toggle Oven on
        """
        self.oven.channel_voltage(self.ovenChannel, 0)
        self.oven.channel_toggle(self.ovenChannel, 0)

    @rpc
    def set_oven_voltage(self, voltage):

        """
        Set the oven voltage
        """
        if not 0. <= voltage <= 5.0:
            raise Exception(f"voltage must be between 0V and 5V")

        self.oven.channel_voltage(self.ovenChannel, voltage)

    @rpc
    def set_oven_current(self, current):

        """
        Set the oven voltage
        """
        if not 0. <= current <= 4.0:
            raise Exception("current must be between 0A and 4A")

        self.oven.channel_current(self.ovenChannel, current)

    @rpc
    def get_oven_voltage(self) -> TFloat:

        return self.oven.measure_voltage(self.ovenChannel)

    @rpc
    def get_oven_current(self) -> TFloat:
        return self.oven.measure_current(self.ovenChannel)
