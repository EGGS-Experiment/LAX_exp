from artiq.experiment import *

import logging
from abc import ABC, abstractmethod
from inspect import getmembers, ismethod

from LAX_exp.base import LAXEnvironment
logger = logging.getLogger("artiq.master.experiments")


class LAXDevice(LAXEnvironment, ABC):
    """
    A generic device object.
        Gets devices and their relevant parameters from the master
        and stores them within the experiment result.

    Attributes:
        name                    str                     : the name of the device (will be used to refer to the device by other objects).
        core_device             tuple(str, str)         : a tuple of (<device_nickname>, <actual_device_name>). All instance methods from
                                                            this device will also become methods of this LAXDevice object.
        devices                 dict(str, str)          : dict of devices used where the key is the device nickname (e.g. "beam_854")
                                                            and the value is the actual device name (e.g. "urukul1_ch3").
    """
    # Class attributes
    name =          None
    core_device =   None
    devices =       dict()


    '''
    BUILD
    '''

    def build(self, **kwargs):
        """
        Get core devices and their parameters from the master, and instantiate them.

        Will be called upon instantiation.
        """
        # get core devices
        self.setattr_device("core")

        # store arguments passed during init for later processing
        self._build_arguments = kwargs

        # set core device as class attribute and break out its methods
        if self.core_device is not None:
            device_nickname, device_name = self.core_device

            try:
                # set device as class attribute
                device_object = self.get_device(device_name)
                setattr(self, device_nickname, device_object)
                self.kernel_invariants.add(device_nickname)

                # verifies that a function is not magic (i.e. a special function, e.g. "__dir__")
                isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (
                            func_obj.__name__ != "__init__")

                # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
                for (function_name, function_object) in getmembers(device_object, isDeviceFunction):
                    setattr(self, function_name, function_object)

            except Exception as e:
                # logger.warning("Device unavailable: {:s}".format(device_name))
                raise Exception("Device unavailable: {:s}".format(device_name))

        # get device(s) & set them as class attributes
        for device_nickname, device_name in self.devices.items():
            # set device as class attribute
            try:
                device_object = self.get_device(device_name)
                setattr(self, device_nickname, device_object)
                self.kernel_invariants.add(device_nickname)
            except Exception as e:
                # logger.warning("Device unavailable: {:s}".format(device_name))
                raise Exception("Device unavailable: {:s}".format(device_name))

        # call user-defined build function
        self.build_device()

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

    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.

        Will be called by parent classes.
        """
        self.prepare_device()

    def prepare_device(self):
        """
        To be subclassed.

        Called after _prepare_parameters.
        Used to customize the device class and set up hardware.
        """
        pass


    '''
    INITIALIZE
    '''

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        """
        To be subclassed.

        todo: document
        """
        pass


    '''
    CLEANUP
    '''

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        """
        To be subclassed.

        todo: document
        """
        pass
