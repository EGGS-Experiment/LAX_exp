from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ

from time import sleep
# todo: move magic/hardcoded configs to device_db_ext


class Flipper(LAXDevice):
    """
    High-level API functions for using the flipper mirror (Thorlabs MFF101)
        to switch between the PMT and the Andor camera.
    Controlled by a Labjack TTL.
    """
    name = "flipper"

    def prepare_device(self):
        # MAGIC NUMBERS
        self.wait_time = 1  # wait time in seconds
        self.time_flipper_trigger = 0.01
        self.flipper_port_name = "EIO7"

        # create labrad connections
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.labjack = self.cxn.labjack_server

        # todo: document
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
            self._flipper_pulse()

        except Exception as e:
            # attempt device reconnection if failure
            if self.reconnect_to_labjack():
                self._flipper_pulse()
            else:
                raise ConnectionError("Cannot Connect to Labjack")

    @rpc
    def _flipper_pulse(self) -> TNone:
        """
        Pulse the flipper via LabJack TTL.
        """
        # ensure TTL is initialized to 0
        self.labjack.write_name(self.flipper_port_name, 0)
        sleep(self.time_flipper_trigger)
        # pulse TTL to flip flipper
        self.labjack.write_name(self.flipper_port_name, 1)
        sleep(self.wait_time)
        self.labjack.write_name(self.flipper_port_name, 0)

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
