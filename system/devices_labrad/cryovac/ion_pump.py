from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class IonPump(LAXDevice):
    """
    High-level API functions for using the ion pump
    """
    name = "ion_pump"
    computer = 'CausewayBay'
    com_port = 'COM9'
    kernel_invariants = {
        "cxn",
    }

    def prepare_device(self):
        self.ion_pump = self.cxn.niops03_server

        device_info = self.ion_pump.device_info()

        computer = device_info[0]
        com_port = device_info[1]

        if computer == '' or com_port == '':
            try:
                self.ion_pump.device_select(self.computer,
                                            self.com_port)
            except Exception as e:
                repr(e)


    @rpc
    def get_status(self) -> TStr:
        return self.ion_pump.status()

    @rpc
    def get_pressure(self) -> TFloat:
        return self.ion_pump.ip_pressure()

    @rpc
    def get_voltage(self) -> TFloat:
        return self.ion_pump.ip_voltage()

    @rpc
    def get_temperature(self) -> TTuple([TFloat, TFloat]):
        return self.ion_pump.ip_temperature()
