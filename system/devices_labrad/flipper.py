import time

from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Flipper(LAXDevice):
    """
    High-level API functions for using the flipper between PMT and Andor Camera
    """
    name = "flipper"

    def prepare_device(self):

        self.wait_time = 1
        self.time_flipper_trigger = 0.01

        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.labjack = self.cxn.labjack_server

        self.flipper_port_name = "DIO3"

        device_handle = self.labjack.device_info()
        if device_handle == -1:
            # get device list
            dev_list = self.labjack.device_list()
            # assume desired labjack is first in list
            self.labjack.device_select(dev_list[0])

    @rpc
    def flip(self) -> TNone:
        """
        Flip the flipper.
        """
        try:
            # ensure labjack is off before sending pulse
            self.labjack.write_name(self.flipper_port_name, 0)
            time.sleep(self.time_flipper_trigger)

            # pulse labjack to flip flipper
            self.labjack.write_name(self.flipper_port_name, 1)
            time.sleep(self.wait_time)
            self.labjack.write_name(self.flipper_port_name, 0)
        except Exception as e:
            if self.reconnect_to_labjack():
                self.labjack.write_name(self.flipper_port_name, 0)
                time.sleep(self.time_flipper_trigger)
                self.labjack.write_name(self.flipper_port_name, 1)
                time.sleep(self.wait_time)
                self.labjack.write_name(self.flipper_port_name, 0)
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
