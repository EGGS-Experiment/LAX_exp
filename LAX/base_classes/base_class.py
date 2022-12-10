from artiq.experiment import *

import logging
from abc import ABC, abstractmethod
from inspect import getmembers, ismethod


logger = logging.getLogger("artiq.master.experiments")


class LAXBaseClass(HasEnvironment, ABC):
    """
    A generic LAX class.
        Sets up an object

    Attributes:
        name                    str                     : the name of the device (will be used to refer to the device by other objects).
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

    # Attributes
    parameters =            dict()


    # SETUP - BUILD
    def build(self):
        """
        Get & instantiate relevant devices from the master and break out their methods.
        """
        # get core device
        self.setattr_device("core")


    # SETUP - PREPARE
    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.
        """
        self._prepare_parameters()
        self.prepare_class()
        self.prepare_hardware()

    def _prepare_parameters(self):
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


    # USER FUNCTIONS
    def prepare_class(self):
        """
        To be subclassed.
        Called after build.
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
