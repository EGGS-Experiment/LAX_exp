"""
Base classes for building experiment sequences for EGGS.
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


    def setattr_argument(self, key, *args, **kwargs):
        # call original setattr_argument function
        super().setattr_argument(key, *args, **kwargs)

        # add key to parameter list to be added to dataset manager later
        self.local_parameters.append(key)


    def prepare(self):
        # set core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # assign global parameters
        for param_name in self.global_parameters:
            # get parameter from dataset manager
            param_value = self._HasEnvironment__dataset_mgr.ddb.get(param_name)
            # set as parameter in dataset manager and HDF5 file
            self._HasEnvironment__dataset_mgr.set(param_name, param_value, archive=False, parameter=True)

        # assign local parameters
        for param_name in self.local_parameters:
            # get parameter from object attributes and store
            param_value = getattr(self, param_name)
            if type(param_value) in (RangeScan, CenterScan):
                param_value = list(param_value)

            self._HasEnvironment__dataset_mgr.set(param_name, param_value, archive=False, parameter=True)

        # call subclassed prepare functions
        self.call_child_method("prepare_experiment")


    @kernel
    def run_reset(self):
        self.core.reset()


    def run(self):
        # reset to remove leftovers from previous experiments
        self.run_reset()

        # prepare devices
        self.call_child_method("prepare_devices")

        # prepare dma
        self.call_child_method("prepare_dma")

        # run experiment
        self.call_child_method("run_experiment")

        # wrap up experiment
        self.finish()


    @kernel
    def finish(self):
        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)
