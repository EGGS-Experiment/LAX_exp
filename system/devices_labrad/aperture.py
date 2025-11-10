from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ
from time import sleep


class Aperture(LAXDevice):
    """
    High-level API functions for using the Aperture Server.
    """
    name = "aperture"
    kernel_invariants = {
        "cxn", "aperture",
        "position_close"
    }

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.aperture = self.cxn.elliptec_server

        self.position_close = 2888  # aperture position for "closed" state

    @rpc
    def open_aperture(self) -> TNone:
        """
        Opens the Aperture
        """
        self.aperture.move_home()

    @rpc
    def close_aperture(self) -> TNone:
        """
        Closes the Aperture
        """
        self.aperture.move_absolute(self.position_close)

    @rpc
    def pulse_aperture_open(self, wait_time_s: TFloat) -> TNone:
        """
        Pulse Aperture Open
        :param wait_time_s: seconds to wait before closing aperture
        """
        self.aperture.move_home()
        sleep(wait_time_s)
        self.aperture.move_absolute(self.position_close)

    @rpc
    def pulse_aperture_close(self, wait_time_s: TFloat) -> TNone:
        """
        Pulse Aperture Close
        :param wait_time_s: seconds to before reopening aperture
        """
        self.aperture.move_absolute(self.position_close)
        sleep(wait_time_s)
        self.aperture.move_home()

