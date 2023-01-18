from artiq.experiment import *

import logging
from numpy import array
from abc import ABC, abstractmethod
from LAX_exp.base.manager_wrappers import LAXDeviceManager, LAXDatasetManager

logger = logging.getLogger("artiq.master.experiments")


class LAXEnvironment(HasEnvironment, ABC):
    """
    Base class for LAX.

    Defines structures and methods common to all core modules in LAX.
    Designed to be used as an LAX version of HasEnvironment.

    Attributes:
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
        kernel_invariants           set(str)                : list of attribute names that won't change while kernel is running.
        devices                     list(LAXDevice)         : list of devices used by the subsequence.
    """
    # Class attributes
    kernel_invariants =         set()


    def __init__(self, managers_or_parent, *args, **kwargs):
        # define instance variables
        self._set_arguments = dict()
        self._build_arguments = dict()

        # correctly extract manager objects
        device_mgr, dataset_mgr, argument_mgr, scheduler_defaults = (None, None, None, None)
        if isinstance(managers_or_parent, tuple):
            device_mgr = managers_or_parent[0]
            dataset_mgr = managers_or_parent[1]
            argument_mgr = managers_or_parent[2]
            scheduler_defaults = managers_or_parent[3]
        # todo: better document; this is for case if we try to use LAX stuff in isolation with regular hasenv
        elif isinstance(managers_or_parent, HasEnvironment):
            device_mgr = managers_or_parent._HasEnvironment__device_mgr
            dataset_mgr = managers_or_parent._HasEnvironment__dataset_mgr
            argument_mgr = managers_or_parent._HasEnvironment__argument_mgr
            scheduler_defaults = {}
            managers_or_parent.register_child(self)
        # todo: better document; this is for case everyone is LAX
        elif isinstance(managers_or_parent, LAXEnvironment):
            device_mgr = managers_or_parent.__device_mgr
            dataset_mgr = managers_or_parent.__dataset_mgr
            argument_mgr = managers_or_parent.__argument_mgr
            scheduler_defaults = {}
            managers_or_parent.register_child(self)

        # wrap manager objects
        if not isinstance(device_mgr, LAXDeviceManager):
            device_mgr = LAXDeviceManager(device_mgr, self)
        if not isinstance(dataset_mgr, LAXDatasetManager):
            dataset_mgr = LAXDatasetManager(dataset_mgr, self)

        # save wrapped manager objects as private variables
        self.__device_mgr = device_mgr
        self.__dataset_mgr = dataset_mgr
        self.__argument_mgr = argument_mgr
        self.__scheduler_defaults = scheduler_defaults

        # initialize parent so we can use their functions which we haven't
        # already subclassed
        managers = (device_mgr, dataset_mgr, argument_mgr, scheduler_defaults)
        super().__init__(managers, *args, **kwargs)


    '''
    ARGUMENTS
    '''

    def setattr_argument(self, *args, **kwargs):
        """
        todo: document
        """
        super().setattr_argument(*args, **kwargs)

        # add argument to _set_arguments (will grab after prepare)
        key, processor = args
        self._set_arguments[key] = None

    def save_arguments(self):
        """
        todo: document
        """
        # add arguments to the dataset manager
        for arg_key in self._set_arguments.keys():
            try:
                # retrieve argument value
                arg_val = getattr(self, arg_key)

                # convert scan objects into an array to ensure compatibility
                if issubclass(type(arg_val), ScanObject):
                    arg_val = array(list(arg_val))

                # store in dataset manager
                self.__dataset_mgr.set(arg_key, arg_val, archive=False, parameter=False, argument=True)
            except KeyError:
                logger.warning("Argument unavailable: {:s}".format(arg_val))


    '''
    PARAMETERS
    '''

    def get_parameter(self, key, group=None, override=True, conversion_function=None, units=None):
        """
        todo: document
        """
        parameter_value = None

        try:
            # get parameter from dataset manager
            if group is not None:
                parameter_key_full = '{}.{}'.format(group, key)
                parameter_value = self.__dataset_mgr.get(parameter_key_full, False, parameter=True)
            else:
                parameter_value = self.__dataset_mgr.get(key, False, parameter=True)

            # get parameter from _build_arguments
            if (override == True) and (key in self._build_arguments):
                parameter_value = self._build_arguments[key]
                self.__dataset_mgr.set(key, parameter_value, archive=False, parameter=True, argument=False)

            # convert parameter as necessary
            if conversion_function is not None:
                if units is None:
                    parameter_value = conversion_function(parameter_value)
                else:
                    parameter_value = conversion_function(parameter_value * units)

            # return parameter
            return parameter_value

        except Exception as e:
            raise Exception("Parameter unavailable: {}".format(key))

    def setattr_parameter(self, key, group=None, override=True, conversion_function=None, units=None):
        """
        # todo: document
        """
        setattr(self, key, self.get_parameter(key, group, override, conversion_function, units))
        self.kernel_invariants.add(key)
