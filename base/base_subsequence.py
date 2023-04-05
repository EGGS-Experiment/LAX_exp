import logging
from numpy import int64, int32
from abc import ABC, abstractmethod

from artiq.experiment import *
from LAX_exp.base import LAXEnvironment
logger = logging.getLogger("artiq.master.experiments")


class LAXSubsequence(LAXEnvironment, ABC):
    """
    Base class for subsequence objects.

    Defines a single short and regularly used pulse sequence to be recorded onto DMA.
    Assumes that the relevant devices have already been initialized.

    Attributes:
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
    """
    # Class attributes
    name = None

    def __init__(self, managers_or_parent, *args, **kwargs):
        # get unique subsequence number to distinguish different DMA handles from the same class
        _dma_count = 0
        if (not isinstance(managers_or_parent, tuple)) and (hasattr(managers_or_parent, '_dma_count')):
            _dma_count = getattr(managers_or_parent, '_dma_count')
            managers_or_parent._dma_count += 1

        # register our dma number for later reference
        setattr(self, '_dma_count', _dma_count)

        # complete regular initialization
        super().__init__(managers_or_parent, *args, **kwargs)


    '''
    BUILD
    '''

    def build(self, **kwargs):
        """
        General construction of the subsequence object.
        Gets/sets instance attributes and process build arguments.

        Called upon instantiation.
        """
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # store arguments passed during init for later processing
        self._build_arguments = kwargs

        # set instance variables
        setattr(self,   'dma_name',                 '{:s}_{:d}'.format(self.name, self._dma_count))
        setattr(self,   'dma_handle',               (0, int64(0), int32(0)))
        setattr(self,   '_dma_record_flag',         False)

        # call user-defined build function
        self.build_subsequence()

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

    def prepare(self):
        """
        Get and convert parameters from the master for use by the device,
        define object methods, and set up the device hardware.

        Will be called by parent classes.
        """
        # store arguments passed during init for later processing
        self._save_arguments()

        # call user-defined prepare function
        self.prepare_subsequence()

    def prepare_subsequence(self):
        """
        To be subclassed.

        Called after _prepare_parameters.
        Used to customize this class.
        """
        pass


    '''
    DMA
    '''

    @kernel(flags={"fast-math"})
    def record_dma(self):
        """
        Records the run subsequence onto core DMA and sets the trace name as an instance attribute.
        """
        # record sequence
        with self.core_dma.record(self.dma_name):
            self.run()

        # set dma record flag
        self._dma_record_flag = True
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def _load_dma(self):
        if self._dma_record_flag == True:
            self.dma_handle = self.core_dma.get_handle(self.dma_name)

    @kernel(flags={"fast-math"})
    def run_dma(self):
        """
        Runs the core sequence from DMA.
        Requires record_dma to have been previously run.
        """
        self.core_dma.playback_handle(self.dma_handle)


    '''
    RUN
    '''

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
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
        Since subsequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass
