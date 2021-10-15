"""
### BEGIN NODE INFO
[info]
name = yzd1
version = 1.0.0
description = Talks to the Lakeshore 336 Temperature Controller
instancename = yzd1

[startup]
cmdline = %PYTHON% %FILE%
timeout = 20

[shutdown]
message = 987654321
timeout = 20
### END NODE INFO
"""

SERVERNAME = 'lakeshore336Server'
TIMEOUT = 1.0
BAUDRATE = 57600
BYTESIZE = 7
STOPBITS = 1
INPUT_CHANNELS = ['A', 'B', 'C', 'D', '0']
OUTPUT_CHANNELS = [1, 2, 3, 4]
TERMINATOR = '\r\n'

from artiq.experiment import *
import numpy as np

from labrad.support import getNodeName
from labrad.devices import LabradServer
from labrad.server import setting

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

class Lakeshore3362(HasEnvironment, LabradServer):
    name = "lakeshore2"
    regKey = "lakeshore2"
    serNode = getNodeName()

    def build(self):
        print("thkimHE")

    # TEMPERATURE DIODES
    @setting(111,'Read Temperature')
    def temperature_read(self, c):
        """
        Get sensor temperature
        """
        print("thkim1")
        #self.setattr_device("ttl7")
        print('thkim2')
        #self.get_device_db()

class TH1(EnvExperiment):
    """TH1"""
    # def run(self):
    #     simple_server_loop({"Lakeshore336": Lakeshore336}, "localhost", 7921)
    def build(self):
        print(self.get_device_db())

    def run(self):
        from labrad import util
        #util.runServer(Lakeshore3362(self))
        util.runServer(yzd1(self))

    def okth1(self):
        print("this works ok")


"""
### BEGIN NODE INFO
[info]
name = yzd1
version = 1.0.0
description = horrible
instancename = yzd1

[startup]
cmdline = %PYTHON% %FILE%
timeout = 20

[shutdown]
message = 987654321
timeout = 20
### END NODE INFO
"""
class yzd1(LabradServer):
    name = "yzd1"
    regKey = "yzd1"
    serNode = getNodeName()

    def __init__(self, ee):
        self.ee = ee
        LabradServer.__init__(self)

    @setting(123, "thkim1")
    def thkim1(self, c):
        "func ok"
        print("hate")
        self.ee.okth1()