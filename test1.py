import numpy as np

from artiq.experiment import *
from sipyco.pc_rpc import simple_server_loop
from sipyco.pc_rpc import Server
from asyncio import get_event_loop

class dummy_device(HasEnvironment):
    def build(self):
        print("thkimHE")

    def func1(self, c):
        """Get sensor temperature"""
        print(c)
        print('thkim2')

class TH1(EnvExperiment):
    """TH1"""

    def build(self):
        print(self.get_device_db())
        self.target = dummy_device(self)
        self.loop = get_event_loop()
        self.server = Server({"device": self.target}, description=None, builtin_terminate=True, allow_parallel=True)

    def run(self):
        self.loop.run_until_complete(self.server.start('::1', 7921))
        try:
            self.loop.run_until_complete(self.server.wait_terminate())
        finally:
            self.loop.run_until_complete(self.server.stop())
        #self.server = simple_server_loop({'device': self.target}, '::1', 7291)

    def okth1(self):
        print("this works ok")