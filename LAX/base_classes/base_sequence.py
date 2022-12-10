from artiq.experiment import *

import logging
from abc import ABC, abstractmethod


logger = logging.getLogger("artiq.master.experiments")


class LAXSequence(HasEnvironment, ABC):
    """
    Base class for sequence objects.
        A set of pulse sequences that achieves a complete ***thing? todo*** (e.g. sideband cooling).
        This should be the only level that compiles onto core DMA.
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
    devices =                       list()
    parameters =                    dict()


    # SETUP
    def build(self):
        """
        Get core devices and their parameters from the master, and instantiate them.
        """
        # get core device
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


    # SETUP - PREPARE
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

    @kernel(flags='fast-math')
    def record_dma(self):
        """
        Record the run sequence onto core DMA.
            This is only offered for convenience, since Sequence objects
            should be the only ones which compile into core DMA.
        Returns:
            str: the handle name.
        """
        # record sequence
        with self.core_dma.record(self.name):
            self.run()

        # return handle
        return self.core_dma.get_handle(self.name)


    # PREPARE - USER FUNCTIONS
    def prepare_class(self):
        """
        To be subclassed.
        Called after _prepare_parameters.
        Used to customize this class.
        """
        pass


    # RUN - USER FUNCTIONS
    @abstractmethod
    def run(self):
        """
        To be subclassed.
        Runs the main pulse sequence.
        Should be a kernel function.
        """
        pass
