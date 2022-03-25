from artiq.experiment import *
from sipyco.pc_rpc import Client
from artiq.master.databases import DatasetDB
from artiq.master.worker_db import DatasetManager
from time import sleep

"""
ds = Client('::1', 3251, 'master_dataset_db')
ds2 = ds.get_dataset('pmt_test_dataset')
print(ds2)
sleep(2)
"""

dataset_db = DatasetDB('C:\\Users\\EGGS1\\Documents\\ARTIQ\\artiq-master\\dataset_db.pyon')
ds = DatasetManager(dataset_db)
ds2 = ds.get('pmt_test_dataset', archive=False)
print(ds2)
sleep(2)