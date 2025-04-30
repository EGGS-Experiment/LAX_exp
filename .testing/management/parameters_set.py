import os
import copy
import h5py
from artiq.experiment import *


class ParametersSet(EnvExperiment):
    """
    Parameters Set
    Used to restore old parameters to the dataset manager in case of problems with dataset_db.pyon.
    """

    def build(self):
        pass

    def run(self):
        file_directory = "C:\\Users\\EGGS1\\Documents\\exp_tmp"

        # change into desired directory and get all h5 files
        os.chdir(file_directory)
        filepath_list = [os.path.join(file_directory, filename)
                         for filename in os.listdir()
                         if ".h5" in filename]


        for filepath in filepath_list:

            # extract values from the given h5 file
            with h5py.File(filepath, 'r') as file:

                # store parameters (these exist only for LAX_exp-type h5 files)
                try:
                    param_dict = copy.deepcopy(dict(file['parameters'].attrs))

                    for param_key, param_val in param_dict.items():
                        try:
                            self.set_dataset(param_key, param_val, broadcast=True, persist=True)
                        except Exception as e:
                            print(e)

                except KeyError:
                    pass

                # store archive parameters (how parameters are normally stored in base ARTIQ)
                try:
                    archive_dict = dict(file['archive'])

                    for param_key, param_val in archive_dict.items():
                        try:
                            self.set_dataset(param_key, param_val[()], broadcast=True, persist=True)
                        except Exception as e:
                            print(e)
                except KeyError:
                    pass


