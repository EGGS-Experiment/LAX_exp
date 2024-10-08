from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Aperture(LAXDevice):
    """
    High-level API functions for using the Aperture Server.
    """
    name = "aperture"

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.aperture = self.cxn.elliptec_server

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
        self.aperture.move_absolute(2888)

    @rpc
    def pulse_aperture_open(self, wait_time: TFloat) -> TNone:
        """
        Pulse Aperture Open
        Args:
            wait_time: seconds to wait before closing aperture again
        """
        self.aperture.move_home()
        time.sleep(wait_time)
        self.aperture.move_absolute(2888)

    @rpc
    def pulse_aperture_close(self, wait_time: TFloat) -> TNone:
        """
        Pulse Aperture Open
        Args:
            wait_time: seconds to before reopening aperture
        """
        self.aperture.move_absolute(2888)
        time.sleep(wait_time)
        self.aperture.move_home()

