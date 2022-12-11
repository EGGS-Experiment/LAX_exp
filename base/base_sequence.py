from artiq.experiment import *

import logging
from abc import ABC, abstractmethod

logger = logging.getLogger("artiq.master.experiments")


class LAXSequence(HasEnvironment, ABC):
    """
    Base class for sequence objects.

    A longer, more complex sequence of operations (compared to Subsequences) and can
        contain non-kernel functions.
    This should not be compiled onto core DMA.
    Assumes that the relevant devices have already been initialized.

    Attributes:
        kernel_invariants           set(str)                : list of attribute names that won't change while kernel is running
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
        devices                     list(LAXDevice)         : list of devices used by the Sequence.
        parameters                  dict(str, (str, func)   : a dict of device parameters. The key will serve as the attribute name,
                                                            and the value is a tuple of the parameter name as stored in dataset_db,
                                                            together with a conversion function (None if no conversion needed).
    """
    # Core attributes
    kernel_invariants =             set()

    # Class attributes
    name =                          None
    parameters =                    dict()
    devices =                       list()


    '''
    BUILD
    '''

    # BUILD - BASE
    def build(self, **kwargs):
        """
        Get core devices and their parameters from the master, and instantiate them.

        Will be called upon instantiation.
        """
        self._build_set_attributes()
        self.build_process_kwargs(**kwargs)
        self._build_process_parameter_kwargs(**kwargs)
        self._build_get_devices()

    def _build_set_attributes(self):
        """
        Set instance attributes.
        """
        setattr(self,   'build_parameters',     dict())

    def _build_process_parameters(self, **kwargs):
        """
        Process kwargs used to initialize the sequence.
        """
        # check kwargs are valid
        if kwargs is not None:

            # get parameter names
            class_parameters = set([parameter_name.split('.')[-1] for parameter_name, _ in self.parameters.values()])
            build_parameters = set([parameter_name for parameter_name, _ in kwargs.keys()])

            # take kwargs passed to the sequence meant to replace parameter values from the dataset manager
            self.build_parameters = {
                parameter_name: kwargs[parameter_name]
                for parameter_name in class_parameters.intersection(build_parameters)
            }

    def _build_get_devices(self):
        """
        Get LAXDevices as well as core devices necessary for the sequence.
        """
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # get LAXDevices
        for device_object in self.devices:

            # set device as class attribute
            try:
                device_name = device_object.name
                setattr(self, device_name, device_object)
                self.kernel_invariants.add(device_name)
            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))


    # BUILD - USER FUNCTIONS
    def build_process_kwargs(self, **kwargs):
        """
        To be subclassed.

        Called after _build_set_attributes.
        Used to process kwargs in a way specific to the sequence.
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
        self._prepare_parameters()
        self.prepare_class()

    def _prepare_parameters(self):
        """
        Get parameters and convert them for use by the sequence.
        """
        # get sequence parameters
        for parameter_name, parameter_attributes in self.parameters.items():
            _parameter_name_dataset, _parameter_conversion_function = parameter_attributes

            try:
                # get parameter from kwargs or dataset manager
                parameter_value = self.build_parameters.get(
                    _parameter_name_dataset,
                    self.get_dataset(_parameter_name_dataset, archive=True)
                )

                # convert parameter as necessary
                if _parameter_conversion_function is not None:
                    parameter_value = _parameter_conversion_function(parameter_value)

                # set parameter as class attribute
                setattr(self, parameter_name, parameter_value)
                self.kernel_invariants.add(parameter_name)

            except Exception as e:
                logger.warning("Parameter unavailable: {:s}".format(_parameter_name_dataset))


    # PREPARE - USER FUNCTIONS
    def prepare_class(self):
        """
        To be subclassed.

        Called after _prepare_parameters.
        Used to customize this class.
        """
        pass


    '''
    RUN
    '''

    # RUN - USER FUNCTIONS
    @abstractmethod
    def run(self):
        """
        Must be subclassed.

        Runs a sequence of operations.
        """
        pass
