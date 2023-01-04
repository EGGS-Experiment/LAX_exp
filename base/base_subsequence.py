from artiq.experiment import *

import logging
from numpy import int64, int32
from abc import ABC, abstractmethod

from LAX_exp.base import LAXBase
logger = logging.getLogger("artiq.master.experiments")


class LAXSubsequence(LAXBase, ABC):
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
    # Class attributes
    devices =                       list()


    def __init__(self, managers_or_parent, *args, **kwargs):
        super().__init__(managers_or_parent, *args, **kwargs)

        # get subseq #
        parent_instance_num_tmp = managers_or_parent.getattr('instance_number', 0)
        setattr(self, 'instance_number', parent_instance_num_tmp)
        print('self inst num: {}'.format(self.instance_number))

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
        self._build_subsequence()
        self.build_subsequence()

    def _build_subsequence(self):
        """
        General construction of the subsequence object.

        Gets/sets instance attributes, devices, and process build arguments.
        Called before build_subsequence.
        """
        # tmp remove
        print('build subseq called')
        # tmp remove clear

        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # set instance variables
        setattr(self,   'dma_name',                 '{:s}_{:d}'.format(self.name, self.instance_number))
        setattr(self,   'dma_handle',               (0, int64(0), int32(0)))
        setattr(self,   '_dma_record_flag',         False)

        # keep track of all instances of a subsequence
        # tmp remove
        print('\tbuild subseq vals set, parent subseq counter: {}'.format(self.instance_number))
        # tmp remove clear

        # set devices as class attributes
        for device_name in self.devices:
            try:
                device_object = self.get_device(device_name)
                setattr(self, device_name, device_object)
                self.kernel_invariants.add(device_name)
            except Exception as e:
                logger.warning("Device unavailable: {:s}".format(device_name))


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
        self._prepare_parameters(**self._build_arguments)
        self.prepare_subsequence()


    # PREPARE - USER FUNCTIONS
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

        Returns:
            str: the DMA handle for the sequence.
        """
        # record sequence
        with self.core_dma.record(self.dma_name):
            self.run()

        # set dma record flag
        self._dma_record_flag = True
        self.core.break_realtime()

    def _load_dma(self):
        if self._dma_record_flag == True:
            self._load_dma_kernel()

            # tmp remove
            print('\t{}: {}'.format(self.name, self.dma_handle))
            print('\t\tdma name: {}'.format(self.dma_name))
            # tmp remove

    # todo: move to solely a kernel function after we fix dma issue
    @kernel(flags={"fast-math"})
    def _load_dma_kernel(self):
        self.dma_handle = self.core_dma.get_handle(self.dma_name)


    '''
    RUN
    '''

    # RUN - BASE
    @kernel(flags={"fast-math"})
    def run_dma(self):
        """
        Runs the core sequence from DMA.
        Requires record_dma to have been previously run.
        """
        self.core_dma.playback_handle(self.dma_handle)
        # todo: remove break_realtime()
        self.core.break_realtime()


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
