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

    # Runtime attributes
    duplicate_counter =             0


    '''
    BUILD
    '''

    # BUILD - BASE
    def build(self, **kwargs):
        """
        Get core devices and their parameters from the master, and instantiate them.

        Will be called upon instantiation.
        """
        self._build_subsequence(**kwargs)
        self.build_subsequence(**kwargs)

    def _build_subsequence(self, **kwargs):
        """
        General construction of the subsequence object.

        Gets/sets instance attributes, devices, and process build arguments.
        Called before build_subsequence.
        """
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # set instance variables
        setattr(self,   'dma_handle',           None)
        setattr(self,   'build_parameters',     dict())
        setattr(self,   'instance_number',      self.duplicate_counter)

        # keep track of all instances of a subsequence
        self.duplicate_counter += 1

        # set devices as class attributes
        for device_name in self.devices:
            try:
                device_object = self.get_device(device_name)
                setattr(self, device_name, device_object)
                self.kernel_invariants.add(device_name)
            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))


        # extract parameters from build arguments
        class_parameters = set([parameter_name.split('.')[-1] for parameter_name, _ in self.parameters.values()])
        build_parameters = set([parameter_name for parameter_name in kwargs.keys()])

        # take kwargs passed to the subsequence meant to replace parameter values from the dataset manager
        self.build_parameters = {
            parameter_name: kwargs[parameter_name]
            for parameter_name in class_parameters.intersection(build_parameters)
        }


    # BUILD - USER FUNCTIONS
    def build_subsequence(self, **kwargs):
        """
        To be subclassed.

        Called after _build_set_attributes.
        Used to set and process subsequence-specific arguments.
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
        self.__prepare_parameters()
        self.prepare_subsequence()
        self._prepare_subsequence()

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
                    self.__dataset_mgr.set(_parameter_name_dataset, parameter_value, archive=False, parameter=True, argument=False)

                # convert parameter as necessary
                if _parameter_conversion_function is not None:
                    parameter_value = _parameter_conversion_function(parameter_value)

                # set parameter as class attribute
                setattr(self, parameter_name, parameter_value)
                self.kernel_invariants.add(parameter_name)

            except Exception as e:
                logger.warning("Parameter unavailable: {:s}".format(_parameter_name_dataset))

    def _prepare_subsequence(self):
        """
        idk
        :return:
        """
        pass


    # PREPARE - USER FUNCTIONS
    def prepare_subsequence(self):
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
    def record_dma(self):
        """
        Records the run subsequence onto core DMA and sets
        the handle as an attribute.
        Returns:
            str: the DMA handle for the sequence.
        """
        # record sequence
        dma_handle = self._record_dma('{:s}_{:d}'.format(self.name, self.instance_number))

        # set dma handle as class attribute
        setattr(self, 'dma_handle', dma_handle)
        self.kernel_invariants.add('dma_handle')

    @kernel(flags='fast-math')
    def _record_dma(self, handle_name):
        # record sequence
        with self.core_dma.record(handle_name):
            self.run()

        # get sequence handle
        self.core.break_realtime()
        handle = self.core_dma.get_handle(handle_name)

        # return handle
        return handle

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
        Since subsequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass


    '''
    HasEnvironment Extensions
    '''

    def get_parameter(self, key, default=NoDefault, archive=False):
        try:
            return self._HasEnvironment__dataset_mgr.get(key, archive, parameter=True)
        except KeyError:
            if default is NoDefault:
                raise
            else:
                return default
