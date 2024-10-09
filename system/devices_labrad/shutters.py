from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Shutters(LAXDevice):
    """
    High-level API functions for using the shutters
    """
    name = "shutters"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.labjack = self.cxn.labjack_server

        self.port_name_377 = "DIO0"
        self.port_name_423 = "DIO2"

        device_handle = self.labjack.device_info()
        if device_handle == -1:
            # get device list
            dev_list = self.labjack.device_list()
            # assume desired labjack is first in list
            self.labjack.device_select(dev_list[0])

    @rpc
    def toggle_377_shutter(self, status: TBool) -> TNone:
        """
        Toggle the 377nm Shutter.
        """
        try:
            self.labjack.write_name(self.port_name_377, status)
        except Exception as e:
            if self.reconnect_to_labjack():
                self.labjack.write_name(self.port_name_377, status)
            else:
                raise ConnectionError("Cannot Connect to Labjack")

    @rpc
    def toggle_423_shutter(self, status: TBool) -> TNone:
        """
        Toggle the 423nm Shutter.
        """
        try:
            self.labjack.write_name(self.port_name_423, status)
        except Exception as e:
            if self.reconnect_to_labjack():
                self.labjack.write_name(self.port_name_423, status)
            else:
                raise ConnectionError("Cannot Connect to Labjack")

    @rpc
    def reconnect_to_labjack(self) -> TNone:
        """
        reconnect to the labjack
        """
        self.labjack.device_close()
        device_handle = self.labjack.device_info()

        if device_handle == -1:
            # get device list
            dev_list = self.labjack.device_list()
            # assume desired labjack is first in list
            self.labjack.device_select(dev_list[0])

            if len(dev_list) >= 1:
                return True
            else:
                return False

        return True
