"""
Results ***

# tmp remove
This module contains the worker process main() function and the glue code
necessary to connect the global artefacts used from experiment code (scheduler,
device database, etc.) to their actual implementation in the parent master
process via IPC.
"""

import h5py
import os.path


def write_results_lax(exp_obj, exp_params):
    rid = exp_params["rid"]

    # get save directories
    save_dir_list = exp_obj.get_dataset('management.dataset_save_locations', [], archive=False)

    # create HDF5 file and save results
    for save_dir in save_dir_list:

        try:
            # create target filepath
            filename = "{:09}-{}.h5".format(rid, exp_obj.name)
            filedir = os.path.join(save_dir, filename)

            # write data
            with h5py.File(filedir, "w") as f:

                # save data from experiment via the dataset manager of the LAXExperiment
                exp_obj._HasEnvironment__dataset_mgr.write_hdf5(f)

                # store experiment parameters in a separate group as attributes
                experiment_group = f.create_group("experiment")
                for k, v in exp_params.items():
                    experiment_group.attrs[k] = v

        # catch any errors
        except Exception as e:
            print("Warning: unable to create and save file in LAX format: {}".format(e))
