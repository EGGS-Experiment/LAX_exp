from sipyco import pyon
import lmdb

old = pyon.load_file("dataset_db.pyon")
new = lmdb.open("dataset_db.mdb", subdir=False, map_size=2**30)
with new.begin(write=True) as txn:
  for key, value in old.items():
    txn.put(key.encode(), pyon.encode((value, {})).encode())
new.close()