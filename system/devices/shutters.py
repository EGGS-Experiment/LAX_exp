from artiq.experiment import *

from LAX_exp.base import LAXDevice

from os import environ
import labrad


class Shutters(LAXDevice):
    """
    High-level api functions for using the Andor Camera
    """

    name = "shutters"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.labjack = self.cxn.labjack_server


        self.port_name_377 = "DIO0"
        self.port_name_423 = "DIO2"

        device_handle = self.labjack.device_info()
        if device_handle is -1:
            # get device list
            dev_list = self.labjack.device_list()
            # assume desired labjack is first in list
            self.labjack.device_select(dev_list[0])

    @rpc
    def open_377_shutter(self) -> TNone:
        """
        Open the Shutter
        """
        self.labjack.write_name(self.port_name_377, 1)

    @rpc
    def close_377_shutter(self) -> TNone:
        """
        Close the Shutter
        """
        self.labjack.write_name(self.port_name_377, 0)

    @rpc
    def open_423_shutter(self) -> TNone:
        """
        Open the Shutter
        """
        self.labjack.write_name(self.port_name_423, 1)

    @rpc
    def close_423_shutter(self) -> TNone:
        """
        Close the Shutter
        """
        self.labjack.write_name(self.port_name_423, 0)

    @rpc
    def close_labjack(self) -> TNone:
        """
        Close labjack
        """
        self.labjack.device_close()
