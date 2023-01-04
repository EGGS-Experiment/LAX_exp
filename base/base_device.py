from artiq.experiment import *

import logging
from abc import ABC, abstractmethod
from inspect import getmembers, ismethod

from LAX_exp.base import LAXBase
logger = logging.getLogger("artiq.master.experiments")


class LAXDevice(LAXBase, ABC):
    """
    A generic device object.
        Gets devices and their relevant parameters from the master
        and stores them within the experiment result.

    Attributes:
        name                    str                     : the name of the device (will be used to refer to the device by other objects).
        kernel_invariants       set(str)                : list of attribute names that won't change while kernel is running
        parameters              dict(str, (str, func)   : a dict of device parameters. The key will serve as the attribute name,
                                                            and the value is a tuple of the parameter name as stored in dataset_db,
                                                            together with a conversion function (None if no conversion needed).
        core_devices            dict(str, str)          : dict of devices used where the key is the device nickname (e.g. "beam_854")
                                                            and the value is the actual device name (e.g. "urukul1_ch3").
    """
    # Class attributes
    core_devices =          dict()


    '''
    BUILD
    '''

    # BUILD - BASE
    def build(self, **kwargs):
        """
        Get core devices and their parameters from the master, and instantiate them.

        Will be called upon instantiation.
        """
        self._build_arguments = kwargs
        self._build_device()
        self.build_device()

    def _build_device(self):
        """
        General instantiation required for the device.

        Sets class attributes and methods.
        """
        # get core device
        self.setattr_device("core")

        # get device(s) & set them as class attributes
        for device_nickname, device_name in self.core_devices.items():

            # set device as class attribute
            try:
                device_object = self.get_device(device_name)
                setattr(self, device_nickname, device_object)
                self.kernel_invariants.add(device_nickname)
            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))

        # if class only uses one device, break out original device methods
        if len(self.core_devices) == 1:

            # verifies that a function is not magic
            isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")

            # get device
            dev_tmp = getattr(self, list(self.core_devices.keys())[0])

            # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
            for (function_name, function_object) in getmembers(dev_tmp, isDeviceFunction):
                setattr(self, function_name, function_object)


    # BUILD - USER FUNCTIONS
    def build_device(self):
        """
        To be subclassed.

        Called at the end of build.
        Used to set & process arguments, define object attributes, etc.
        """
        pass


    '''
    PREPARE
    '''

    # PREPARE - BASE
    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.

        Will be called by parent classes.
        """
        self._prepare_parameters(**self._build_arguments)
        self.prepare_device()

    # PREPARE - USER FUNCTIONS
    def prepare_device(self):
        """
        To be subclassed.

        Called after _prepare_parameters.
        Used to customize the device class and set up hardware.
        """
        pass
