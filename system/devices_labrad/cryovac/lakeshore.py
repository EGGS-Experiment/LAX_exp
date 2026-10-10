from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ


class Lakeshore(LAXDevice):
    """
    High-level API functions for using the lakeshore
    """
    name = "lakeshore"
    kernel_invariants = {
        "cxn",
    }

    def prepare_device(self):
        self.lakeshore = self.cxn.lakeshore336_server

        device_info = self.lakeshore.device_info()

        computer = device_info[0]
        com_port = device_info[1]

        if computer == '' or com_port == '':
            try:
                self.lakeshore.device_select()
            except Exception as e:
                repr(e)

    def get_temperature(self, diode = "All") -> TList(TFloat):

        if diode not in ['All', 'A', 'B', 'C', 'D']:
            raise ValueError('temperature diode must given as All, A, or B, or C, or D')

        if diode == "All":
            temp_a = self.lakeshore.read_temperature('A')[0]
            temp_b = self.lakeshore.read_temperature('B')[0]
            temp_c = self.lakeshore.read_temperature('C')[0]
            temp_d = self.lakeshore.read_temperature('D')[0]
            temps = [temp_a, temp_b, temp_c, temp_d]

        else:
            temps = [self.lakeshore.read_temperature(diode)]

        return temps




