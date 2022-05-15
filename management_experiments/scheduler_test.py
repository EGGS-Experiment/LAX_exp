from artiq.experiment import *
from sipyco.pc_rpc import Client
from artiq.master.databases import DatasetDB
from artiq.master.worker_db import DatasetManager
from time import sleep

sh1 = Client('::1', 3251, 'master_schedule')
status = sh1.get_status()
print(status)
pause = sh1.check_pause(1000)
print(pause)