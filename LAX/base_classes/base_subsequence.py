from artiq.experiment import *

import logging
from abc import ABC, abstractmethod

logger = logging.getLogger("artiq.master.experiments")


class LAXSubsequence(HasEnvironment, ABC):
    """
    Base class for subsequence objects.
        Defines a single short and regularly used pulse sequence to be recorded onto DMA.
        Assumes that the relevant devices have already been initialized.

    Attributes:
        kernel_invariants           set(str)                : list of attribute names that won't change while kernel is running
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
        devices                     list(LAXDevice)         : list of devices used by the subsequence.
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
        self._build_process_kwargs(**kwargs)
        self._build_get_devices()

    def _build_set_attributes(self):
        """
        Set instance attributes.
        """
        # set instance variables
        setattr(self,   'dma_handle',           None)
        setattr(self,   'build_parameters',     dict())

    def _build_process_kwargs(self, **kwargs):
        """
        Process kwargs used to initialize the subsequence.
        """
        # check kwargs are valid
        if kwargs is not None:

            # get parameter names
            valid_parameters = set([parameter_name for parameter_name, _ in self.parameters.values()])
            build_parameters = set([parameter_name for parameter_name, _ in kwargs.keys()])

            # make passed args
            if build_parameters.issubset(valid_parameters):
                self.build_parameters = kwargs

    def _build_get_devices(self):
        """
        Get LAXDevices as well as core devices necessary for the subsequence.
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
        self._prepare_subsequence()

    def _prepare_parameters(self):
        """
        Get parameters and convert them for use by the subsequence.
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

                # convert parameter to machine units as necessary
                if _parameter_conversion_function is not None:
                    parameter_value = _parameter_conversion_function(parameter_value)

                # set parameter as class attribute
                setattr(self, parameter_name, parameter_value)
                self.kernel_invariants.add(parameter_name)

            except Exception as e:
                logger.warning("Parameter unavailable: {:s}".format(parameter_name))

    def _prepare_subsequence(self):
        """
        Records the subsequence onto core DMA and sets
        the handle as an attribute.
        """
        # record sequence
        dma_handle = self._prepare_dma()

        # set dma handle as class attribute
        setattr(self, 'dma_handle', dma_handle)
        self.kernel_invariants.add(dma_handle)

    @kernel(flags='fast-math')
    def _prepare_dma(self):
        """
        Record the run sequence onto core DMA.
        Returns:
            str: the DMA handle for the sequence.
        """
        # record sequence
        with self.core_dma.record(self.name):
            self.run()

        # get sequence handle
        self.core.break_realtime()
        handle = self.core_dma.get_handle(self.name)

        # return handle
        return handle


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

    # RUN - BASE
    @kernel(flags='fast-math')
    def run_dma(self):
        """
        Runs the core sequence from DMA.
        Requires _prepare_subsequence to have already been run.
        """
        self.core_dma.playback_handle(self.dma_handle)


    # RUN - USER FUNCTIONS
    @abstractmethod
    def run(self):
        """
        Must be subclassed.

        Runs a fixed, unchangeable pulse sequence from core DMA.
        Must have the kernel decorator.
        Since Subsequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass
