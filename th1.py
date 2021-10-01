from artiq.experiment import *
import numpy as np
import time

from serial import Serial, PARITY_ODD
from sipyco.pc_rpc import simple_server_loop
import sipyco.common_args as sca

class Lakeshore336(object):
    def __init__(self, device_name):
        self.ser = Serial(device_name, baudrate = 57600, bytesize = 7, parity = PARITY_ODD, timeout = 1.0)

    def identify(self):
        self.ser.write("*IDN?\n".encode())
        return self.ser.readline().decode()

    def get_temp(self, input='0'):
        """
        Returns the temperature of an input channel as a float in Kelvin
        Args:
            input (str): either 0, 'A', 'B', 'C', or 'D' (0 is all channels)
        """
        self.ser.write('KRDG? ' + str(output_channel) + TERMINATOR)
        resp = self.ser.read_until()
        resp = np.array(resp.split(','), dtype=float)
        print(resp)
        return resp

    def close(self):
        self.ser.close()

class TH1(EnvExperiment):
    """Management Tutorial"""
    def run(self):
        try:
            simple_server_loop({"Lakeshore336": Lakeshore336}, "::1", 3259)
        finally:
            dev.close()