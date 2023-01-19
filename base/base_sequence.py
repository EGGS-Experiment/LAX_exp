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
    # Class attributes
    name = None

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

        # call user-defined build function
        self.build_sequence()

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

    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.

        Will be called by parent classes.
        """
        # store arguments in dataset manager
        self._save_arguments()

        # call user-defined prepare function
        self.prepare_sequence()

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

    @kernel(flags={"fast-math"})
    def initialize_sequence(self):
        """
        To be subclassed.

        todo: document
        note: don't initialize devices here, otherwise lots of redundancy and overhead
        """
        pass

    @abstractmethod
    def run(self):
        """
        Must be subclassed.

        Runs a fixed, unchangeable pulse sequence from core DMA.
        Must have the kernel decorator.
        Since sequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass
