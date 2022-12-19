from artiq.experiment import *

import logging
from abc import ABC, abstractmethod

logger = logging.getLogger("artiq.master.experiments")


class LAXBase(HasEnvironment, ABC):
    """
    Base class for LAX.

    Defines structures and methods common to all core modules in LAX.
    Designed to be used as an LAX version of HasEnvironment..

    Attributes:
        kernel_invariants           set(str)                : list of attribute names that won't change while kernel is running
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
        devices                     list(LAXDevice)         : list of devices used by the subsequence.
        parameters                  dict(str, (str, func)   : a dict of device parameters. The key will serve as the attribute name,
                                                            and the value is a tuple of the parameter name as stored in dataset_db,
                                                            together with a conversion function (None if no conversion needed).
    """
    # Core attributes
    kernel_invariants = set()

    # Class attributes
    name = None


    def __prepare_parameters(self):
        """
        Get parameters and convert them for use by the subsequence.
        """
        # get sequence parameters
        for parameter_name, parameter_attributes in self.parameters.items():
            _parameter_name_dataset, _parameter_conversion_function = parameter_attributes

            try:
                # get parameter from kwargs
                parameter_value = self.build_parameters.get(
                    _parameter_name_dataset,
                    None
                )

                # if not in kwargs, take it from dataset manager
                if parameter_value is None:
                    parameter_value = self.get_parameter(_parameter_name_dataset)

                # if in kwargs, add parameter to dataset manager
                else:
                    self.__dataset_mgr.set(_parameter_name_dataset, parameter_value, archive=False, parameter=True,
                                           argument=False)

                # convert parameter as necessary
                if _parameter_conversion_function is not None:
                    parameter_value = _parameter_conversion_function(parameter_value)

                # set parameter as class attribute
                setattr(self, parameter_name, parameter_value)
                self.kernel_invariants.add(parameter_name)

            except Exception as e:
                logger.warning("Parameter unavailable: {:s}".format(_parameter_name_dataset))

    def get_parameter(self, key, default=NoDefault, archive=False):
        try:
            return self._HasEnvironment__dataset_mgr.get(key, archive, parameter=True)
        except KeyError:
            if default is NoDefault:
                raise
            else:
                return default

    def setattr_argument(self, *args, **kwargs):
        super().setattr_argument(*args, **kwargs)

        # add argument to _build_arguments (will grab after prepare)
        key, processor = args
        self._build_arguments[key] = None
