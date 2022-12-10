from artiq.experiment import *

import logging
from abc import ABC, abstractmethod
from inspect import getmembers, ismethod

logger = logging.getLogger("artiq.master.experiments")


class LAXDevice(HasEnvironment, ABC):
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
    # Core attributes
    name =                  None
    kernel_invariants =     set()

    # Class attributes
    parameters =            dict()
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
        #self._build_set_attributes()
        #self._build_process_kwargs(**kwargs)
        self._build_get_devices()
        self._build_expose_methods()

    def _build_get_devices(self):
        """
        Get core devices necessary for the subsequence.
        """
        # get core device
        self.setattr_device("core")

        # get device(s)
        for device_nickname, device_name in self.core_devices.items():

            # set device as class attribute
            try:
                device_object = self.get_device(device_name)
                setattr(self, device_nickname, device_object)
                self.kernel_invariants.add(device_nickname)
            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))

    def _build_expose_methods(self):
        """
        Expose core device methods as the object's own,
        i.e. allow users to use underlying core device methods directly from this class.
        """
        # if class only uses one device, break out original device methods
        if len(self.core_devices) == 1:

            # verifies that a function is not magic
            isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")

            # get device
            dev_tmp = getattr(self, list(self.core_devices.keys())[0])

            # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
            for (function_name, function_object) in getmembers(dev_tmp, isDeviceFunction):
                setattr(self, function_name, function_object)


    '''
    PREPARE
    '''

    # PREPARE - BASE
    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.

        To be called by parent classes.
        """
        self._prepare_parameters()
        self.prepare_class()
        self.prepare_hardware()

    def _prepare_parameters(self):
        """
        Get parameters and convert them for use by the object.
        """
        # get device parameters
        for parameter_name, parameter_attributes in self.parameters.items():
            _parameter_name_dataset, _parameter_conversion_function = parameter_attributes

            # set parameter as class attribute
            try:
                # get parameter from dataset manager and store in HDF5
                parameter_value = self.get_dataset(_parameter_name_dataset, archive=True)

                # convert parameter to machine units as necessary
                if _parameter_conversion_function is not None:
                    parameter_value = _parameter_conversion_function(parameter_value)

                setattr(self, parameter_name, parameter_value)
                self.kernel_invariants.add(parameter_name)

            except Exception as e:
                logger.warning("Parameter unavailable: {:s}".format(parameter_name))


    # PREPARE - USER FUNCTIONS
    def prepare_class(self):
        """
        To be subclassed.

        Called after _prepare_parameters.
        Used to customize this class.
        """
        pass

    def prepare_hardware(self):
        """
        To be subclassed.

        Called after prepare_class.
        Used to set up device hardware.
        """
        pass
