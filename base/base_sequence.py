import logging
from numpy import int64, int32
from abc import ABC, abstractmethod

from artiq.experiment import *
from LAX_exp.base import LAXEnvironment
logger = logging.getLogger("artiq.master.experiments")


class LAXSequence(LAXEnvironment, ABC):
    """
    Base class for Sequence objects.

    Defines a single short and regularly used pulse sequence.
    Assumes that the relevant devices have already been initialized.

    Attributes:
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
    """

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
        self._build_sequence()
        self.build_sequence()

    def _build_sequence(self):
        """
        General construction of the sequence object.

        Gets/sets instance attributes, devices, and process build arguments.
        Called before build_sequence.
        """
        # get core devices
        self.setattr_device("core")

    # BUILD - USER FUNCTIONS
    def build_sequence(self, **kwargs):
        """
        To be subclassed.

        Called after _build_set_attributes.
        Used to set and process sequence-specific arguments.
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
        self.save_arguments()
        self.prepare_sequence()

    # PREPARE - USER FUNCTIONS
    def prepare_sequence(self):
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

        Runs a fixed, unchangeable pulse sequence from core DMA.
        Must have the kernel decorator.
        Since sequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass
