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
        kernel_invariants       set(str)                : list of attribute names that won't change while kernel is running
        device_names            dict(str, str)          : dict of devices used where the key is the device nickname (e.g. "beam_854")
                                                            and the value is the actual device name (e.g. "urukul1_ch3").
        device_parameters       dict(str, (str, func)   : a dict of device parameters. The key will serve as the attribute name,
                                                            and the value is a tuple of the parameter name as stored in dataset_db,
                                                            together with a conversion function (None if no conversion needed).
    """
    # Core attributes
    name =                  None
    kernel_invariants =     set()

    # Device attributes
    device_names =          dict()
    device_parameters =     dict()

    # todo: document
    _initialized = False


    # SETUP
    def build(self):
        """
        Get core devices and their parameters from the master, and instantiate them.
        """
        # ensure device is only automatically built once per experiment
        try:
            # get list of all extant devices
            _laxDevice_list_all = self.get_dataset('_laxDevice_list_all', set(), archive=False)

            # stop here if device already exists
            if self.name in _laxDevice_list_all:
                return
            # otherwise, add this device to the list
            else:
                _laxDevice_list_all.add(self.name)
                self.set_dataset('_laxDevice_list_all', _laxDevice_list_all, archive=True)

        # begin list of extant devices if it doesn't yet exist
        except KeyError:
            self.set_dataset('_laxDevice_list_all', set([self.name]), archive=True)


        # get core device
        self.setattr_device("core")

        # get device(s)
        for device_nickname, device_name in self.device_names.items():

            # set device as class attribute
            try:
                device_object = self.get_device(device_name)

                setattr(self, device_nickname, device_object)
                self.kernel_invariants.add(device_nickname)

            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))

        # get device parameters
        for parameter_name, parameter_attributes in self.device_parameters.items():

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

        # if class only uses one device, break out original device methods
        if len(self.device_names) == 1:

            # verifies that a function is not magic
            isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")

            # get device
            dev_tmp = getattr(self, list(self.device_names.keys())[0])

            # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
            for (function_name, function_object) in getmembers(dev_tmp, isDeviceFunction):
                setattr(self, function_name, function_object)

        # run device-specific initialization
        self.prepare_devices()


    @abstractmethod
    def prepare_devices(self):
        """
        To be subclassed.
        Called after initClient.
        Used to instantiate a GUI class with arbitrary configuration settings.
        """
        pass
