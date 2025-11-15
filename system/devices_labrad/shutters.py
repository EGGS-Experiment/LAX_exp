from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ

# todo: move magic/hardcoded configs to device_db_ext


class Shutters(LAXDevice):
    """
    High-level API functions for using the shutters
    """
    name = "shutters"
    kernel_invariants = {
        "cxn", "labjack",
        "port_name_377", "port_name_423",
    }

    def prepare_device(self):
        # establish labrad client connections
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.labjack = self.cxn.labjack_server

        # MAGIC NUMBERS
        self.port_name_377 = "DIO0"
        self.port_name_423 = "DIO2"

        # todo: document
        # maybe: use reconnect via labjack?
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
        :param status: 423nm shutter status. True is OPEN, False is OFF.
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
        :param status: 423nm shutter status. True is OPEN, False is OFF.
        """
        try:
            self.labjack.write_name(self.port_name_423, status)
        except Exception as e:
            if self.reconnect_to_labjack():
                self.labjack.write_name(self.port_name_423, status)
            else:
                raise ConnectionError("Cannot Connect to Labjack")

    @rpc
    def reconnect_to_labjack(self) -> TBool:
        """
        Attempt to reconnect to the labjack device.
        :return: device reconnection success status.
        """
        # ensure any existing device connections are closed
        self.labjack.device_close()
        if self.labjack.device_info() == -1:
            # get list of available devices
            dev_list = self.labjack.device_list()

            if len(dev_list) >= 1:
                # assume desired labjack is first in list
                self.labjack.device_select(dev_list[0])
            else:
                return False
        return True
