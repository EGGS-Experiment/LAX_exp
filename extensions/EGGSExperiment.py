"""
Base classes for building PyQt5 GUI clients for LabRAD.
"""
import numpy as np
from artiq.experiment import *

__all__ = ["EGGSExperiment"]


class EGGSExperiment(EnvExperiment):
    """
    A typical EGGS Experiment.
    Does all necessary startup steps under the hood such that the user only
        needs to do implementation-specific initialization.

    Attributes:
        name            (str)               : the name of the GUI client.
        servers         (dict of str: str)  : specifies the LabRAD servers that the client needs and
                                                sets them as instance variables with the name as the key.
        gui             (QWidget)           : the gui object associated with the client.
        LABRADHOST      (str)               : the IP address of the LabRAD manager. If left as None,
                                                it will be set to the LABRADHOST environment variable.
        LABRADPASSWORD  (str)               : the password used to connect to the LabRAD manager. If left as None,
                                                it will be set to the LABRADHOST environment variable.
    """

    # parameter management
    local_parameters = []
    global_parameters = []


    def prepare(self):
        # assign global parameters
        for param_name in self.global_parameters:
            # get parameter from dataset manager
            param_value = self._HasEnvironment__dataset_mgr.ddb.get(param_name)
            # set as parameter in dataset manager and HDF5 file
            self._HasEnvironment__dataset_mgr.set(param_name, param_value, archive=False, parameter=True)
        # assign local parameters
        for param_name in self.local_parameters:
            # get parameter from object attributes
            param_value = getattr(self, param_name)
            # store parameter
            if type(param_value) in (RangeScan, CenterScan):
                param_value = list(param_value)

            self._HasEnvironment__dataset_mgr.set(param_name, param_value, archive=False, parameter=True)

        # call subclassed prepare function
        self.call_child_method("prepare2")

    def setattr_argument(self, key, *args, **kwargs):
        # call original setattr_argument function
        super().setattr_argument(key, *args, **kwargs)

        # add key to parameter list to be added to dataset manager later
        self.local_parameters.append(key)

    # todo: upon removal stuff
