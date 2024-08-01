from artiq.experiment import *

from LAX_exp.base import LAXDevice

from os import environ
import labrad


class Shutters(LAXDevice):
    """
    High-level api functions for using the Andor Camera
    """

    name = "Shutters"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.labjack = self.cxn.labjack_server

        self.port_name_377 = "DIO0"
        self.port_name_423 = "DIO2"

    @kernel(flags={"fast-math"})
    def initialize_device(self):

        if self.labjack.device_handle is not None:
            self.labjack.device_select(["T7", "ANY", "ANY"])

    @rpc
    def open_377_shutter(self, port) -> TNone:
        """
        Open the Shutter
        """
        self.labjack.write_name(self.port_name_377, 1)

    @rpc
    def close_377_shutter(self, port) -> TNone:
        """
        Close the Shutter
        """
        self.labjack.write_name(self.port_name_377, 0)

    @rpc
    def open_423_shutter(self, port) -> TNone:
        """
        Open the Shutter
        """
        self.labjack.write_name(self.port_name_423, 1)

    @rpc
    def close_423_shutter(self, port) -> TNone:
        """
        Close the Shutter
        """
        self.labjack.write_name(self.port_name_423, 0)

    @rpc
    def close_labjack(self) -> TNone:

        self.labjack.device_close()
