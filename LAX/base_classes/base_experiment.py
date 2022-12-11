from artiq.experiment import *

import logging
from abc import ABC, abstractmethod

from LAX_exp.LAX.extensions import LAXDeviceManager, LAXDatasetManager
logger = logging.getLogger("artiq.master.experiments")


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


    '''
    BUILD
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # wrap manager objects
        self.__device_mgr = LAXDeviceManager(self.__device_mgr)
        self.__dataset_mgr = LAXDatasetManager(self.__dataset_mgr)


    # BUILD - BASE
    def build(self):
        """
        Get core devices and their parameters from the master, and instantiate them.

        Will be called upon instantiation.
        """
        self._build_experiment()
        self.build_experiment()

    def _build_experiment(self):
        """
        General construction of the experiment object.

        Gets/sets instance attributes, devices, and process build arguments.
        Called before build_experiment.
        """
        # core devices, etc.
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_cpld')

        # instance variables
        setattr(self,   'dma_handle',               None)
        setattr(self,   '_build_arguments',         dict())


        # universal arguments
        self.setattr_argument("repetitions",                    NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))


    # BUILD - USER FUNCTIONS
    def build_experiment(self):
        """
        To be subclassed.

        Called during build.
        Used to specify experiment arguments and get relevant devices/objects.
        """
        pass


    '''
    PREPARE
    '''

    # PREPARE - BASE
    def prepare(self):
        """
        General construction of the experiment object.

        Called right before run (unlike build, which is called upon instantiation).
        _prepare_experiment is called after prepare_experiment since subsequence instantiation may require
            values computed only in prepare_experiment.
        """
        self.prepare_experiment()
        self._prepare_experiment()

    def _prepare_experiment(self):
        """
        General construction of the experiment object.
        """
        # add arguments to the dataset manager
        for arg_key in self._build_arguments.keys():

            try:
                arg_val = getattr(self, arg_key)
                self.__dataset_mgr.set(arg_key, arg_val, archive=False, parameter=False, argument=True)
            except KeyError:
                logger.warning("Argument unavailable: {:s}".format(arg_val))

        # create dataset to hold results
        self.set_dataset('results', list())
        self.setattr_dataset('results')

        # prepare all children
        self.call_child_method("prepare")


    # PREPARE - USER FUNCTIONS
    def prepare_experiment(self):
        """
        To be subclassed.

        Used to pre-compute data, modify arguments, and further instantiate objects.
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

    @rpc(flags='async')
    def update_dataset(self):
        """
        To be subclassed.

        Updates the experiment dataset.
        Used to move as much processing off the kernel and onto the host as possible.
        Should use the @rpc decorator with the "async" flag for efficiency.
        """
        pass


    '''
     HasEnvironment Extensions
     '''

    def setattr_argument(self, key, **kwargs):
        super().setattr_argument(key, **kwargs)

        # add argument to _build_arguments (will grab after prepare)
        self._build_arguments[key] = None
