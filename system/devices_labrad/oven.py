from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Oven(LAXDevice):
    """
    High-level API functions for using the Oven
    """
    name = "oven"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.oven = self.cxn.gpp3060_server
        self.ovenChannel = 1

    @rpc
    def toggle(self, status: TBool) -> TNone:
        """
        Toggle Oven on/off.
        Args:
            status (TBool): on/off status of oven.
        """
        self.oven.channel_toggle(self.ovenChannel, status)

    @rpc
    def set_oven_voltage(self, voltage: TFloat) -> TNone:
        """
        Set the oven voltage.
        """
        if not 0. <= voltage <= 5.0:
            raise Exception(f"voltage must be between 0V and 5V")

        self.oven.channel_voltage(self.ovenChannel, voltage)

    @rpc
    def set_oven_current(self, current: TFloat) -> TNone:
        """
        Set the oven voltage.
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
