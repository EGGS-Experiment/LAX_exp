from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Flowmeter(LAXDevice):
    """
    High-level API functions for reading the flowmeter
    """
    name = "flowmeter"
    port_name_flowmeter = "AIN0"
    kernel_invariants = {
        "cxn"
    }

    def prepare_device(self):
        self.labjack = self.cxn.labjack_server

        device_handle = self.labjack.device_info()
        if device_handle == -1:
            # get device list
            dev_list = self.labjack.device_list()
            # assume desired labjack is first in list
            self.labjack.device_select(dev_list[0])

        self.labjack.write_name((
            f"{self.port_name_flowmeter}_RANGE",
            10,
        ))

        self.labjack.write_name((
            f"{self.port_name_flowmeter}_RESOLUTION_INDEX",
            12,
        ))

        self.voltage_max_v = 5.
        self.flow_max_lpm = 20.


    @rpc
    def get_flow_lpm(self) -> TFloat:

        try:
           voltage_read =  self.labjack.read_name(self.port_name_flowmeter)
        except Exception as e:
            if self.reconnect_to_labjack():
                voltage_read = self.labjack.read_name(self.port_name_flowmeter)
            else:
                raise ConnectionError("Cannot Connect to Labjack")

        return voltage_read / self.voltage_max_v * self.flow_max_lpm


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

