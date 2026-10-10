from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class F70Compressor(LAXDevice):
    """
    High-level API functions for using the F70 compressor
    """
    name = "f70_compressor"
    computer = 'CausewayBay'
    kernel_invariants = {
        "cxn",
    }

    def prepare_device(self):
        # establish labrad connection
        self.f70_compressor = self.cxn.f70_server

    @rpc
    def turn_on(self) -> TStr:
        return self.f70_compressor.turn_on()

    @rpc
    def turn_off(self) -> TStr:
        return self.f70_compressor.turn_off()

    @rpc
    def get_status(self) -> TStr:
        return self.f70_compressor.status()


