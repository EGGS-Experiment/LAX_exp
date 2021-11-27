from artiq.experiment import *
from artiq.master.databases import DeviceDB
from artiq.master.worker_db import DeviceManager
from sipyco.pc_rpc import Client

from labrad.server import LabradServer, setting
from twisted.internet.threads import deferToThread
# print('th1')
# from artiq.coredevice.ad9910 import AD9910
# ad9910=AD9910()
# print(ad9910.frequency_to_ftw(10000))
print('th2')
devices=DeviceDB('C:\\Users\\EGGS1\\Documents\\ARTIQ\\artiq-master\\device_db.py')
dm=DeviceManager(devices)
core=dm.get('core')
ttl4=dm.get('ttl4')
print(core)
print(ttl4)

class api(object):
    def __init__(self):
        self.core=core

    @kernel
    def on(self):
        core.reset()
        ttl4.on()
        print('thkim')

api_obj=api()
api_obj.on()