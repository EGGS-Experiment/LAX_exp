from artiq.experiment import *

import logging
from abc import ABC, abstractmethod

logger = logging.getLogger("artiq.master.experiments")

# todo: instantiate one of each device
# todo: make device references within seq/subseq use our device object
# todo: modify device manager

# todo: prepare sequences
# todo: prepare devices

# todo:

# todo: make device parameters be stored separately (in a group)
# todo: make sequence parameters be stored separately (in a group)


class LAXExperiment(EnvExperiment, ABC):
    """
    Base class for experiment objects.

    Runs a sequence a given number of times and records corresponding data.
    Instantiates and initalizes all relevant devices, sequences, and subsequences.

    Attributes:
        kernel_invariants           set(str)                : list of attribute names that won't change while kernel is running
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
        devices                     list(LAXDevice)         : list of devices used by the subsequence.
    """
    # Core attributes
    kernel_invariants =             set()

    # Class attributes
    name =                          None
    devices =                       set()


    '''
    BUILD
    '''

    # BUILD - BASE
    def build(self):
        """
        Specify & get arguments, and get & instantiate core devices and necessary modules
        """
        # tmp remove
        self.setattr_argument("repetitions",                    NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.urukul0 = self.get_device('urukul0_cpld')
        self.urukul1 = self.get_device('urukul1_cpld')

        # tmp remove
        self.build_arguments()
        self._build_set_attributes()
        self.build_children()
        self._build_children()

    def _build_set_attributes(self):
        """
        Set instance attributes.
        """
        setattr(self,   'dma_handle',           None)
        setattr(self,   'build_parameters',     dict())

    def _build_children(self):
        """
        todo: document properly
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


    # BUILD - USER FUNCTIONS
    def build_arguments(self):
        """
        To be subclassed.

        Called during build.
        Used to specify experiment arguments using setattr_argument.
        """
        pass

    def build_children(self):
        """
        To be subclassed.

        Called during _build_set_attributes.
        Used to get relevant core devices, LAX devices, sequences, and subsequences.
        """
        pass


    '''
    PREPARE
    '''

    # PREPARE - BASE
    def prepare(self):
        """
        Get and convert arguments, and set up child devices/sequences.
        """
        self.prepare_arguments()
        self.prepare_children()
        self._prepare_children()
        self._prepare_datasets()

    def _prepare_children(self):
        """
        Collate all required children and run child-specific setup.
        """
        pass

    def _prepare_datasets(self):
        """
        Prepare datasets for use.
        """
        pass


    # PREPARE - USER FUNCTIONS
    def prepare_arguments(self):
        """
        To be subclassed.

        Called after prepare.
        Used to modify/prepare experiment arguments from build.
        """
        pass

    def prepare_children(self):
        """
        To be subclassed.

        Called after prepare_arguments.
        Used to declare and instantiate child objects.
        """
        pass


    '''
    RUN
    '''

    # RUN - BASE
    def run(self):
        """
        Main sequence of the experiment.
        Repeat a given sequence a number of times.
        """
        try:
            # repeat the experimental sequence
            for trial_num in range(self.repetitions):

                # prepare the trial
                self.run_prepare()

                # run the trial
                self.run_loop()

        except TerminationRequested:
            pass

        finally:
            # set devices back to their default state
            self._run_cleanup()

    @kernel(flags='fast-math')
    def _run_cleanup(self):
        """
        Set all devices back to their original state.
        """
        self.urukul0.set_profile(0)
        self.urukul1.set_profile(0)


    # RUN - USER FUNCTIONS
    def run_prepare(self):
        """
        To be subclassed.

        Runs a fixed, unchangeable pulse sequence from core DMA.
        Must have the kernel decorator.
        Since subsequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass

    def run_loop(self):
        """
        Must be subclassed.

        Runs a fixed, unchangeable pulse sequence from core DMA.
        Must have the kernel decorator.
        Since subsequences use core DMA, it cannot contain any methods involving RTIO input.
        """
        pass

    def update_dataset(self):
        """
        To be subclassed.

        Updates the experiment dataset.
        Used to move as much processing off the kernel and onto the host as possible.
        Should use the @rpc decorator with the "async" flag for efficiency.
        """
        pass
