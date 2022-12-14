from artiq.experiment import *

import logging
from numpy import array
from abc import ABC, abstractmethod
logger = logging.getLogger("artiq.master.experiments")

from LAX_exp.extensions import LAXDeviceManager, LAXDatasetManager


class LAXExperiment(HasEnvironment, ABC):
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
        self._HasEnvironment__device_mgr = LAXDeviceManager(self._HasEnvironment__device_mgr, self)
        self._HasEnvironment__dataset_mgr = LAXDatasetManager(self._HasEnvironment__dataset_mgr, self)

        # tmp remove
        self.__device_mgr = self._HasEnvironment__device_mgr
        self.__dataset_mgr = self._HasEnvironment__dataset_mgr

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
        setattr(self,   '_build_arguments',         dict())
        setattr(self,   '_result_iter',             0)

        # universal arguments
        self.setattr_argument("repetitions",                    NumberValue(default=4, ndecimals=0, step=1, min=1, max=10000))


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
        Must happen after the user-defined prepare_experiment method.
        """
        # add arguments to the dataset manager
        for arg_key in self._build_arguments.keys():

            try:
                arg_val = getattr(self, arg_key)
                self.__dataset_mgr.set(arg_key, arg_val, archive=False, parameter=False, argument=True)
            except KeyError:
                logger.warning("Argument unavailable: {:s}".format(arg_val))

        # create dataset to hold results
        #self.set_dataset('results', list())
        #self.setattr_dataset('results')

        # prepare children
        self.call_child_method('prepare')


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
        # set up the run
        self.run_initialize()

        try:
            # repeat the experiment a given number of times
            for trial_num in range(self.repetitions):

                # prepare a trial
                self.loop_prepare()

                # run the trial
                self.loop_run()

        # allow clean termination
        except TerminationRequested:
            pass

        # set devices back to their default state
        finally:
            self._run_cleanup()

    @kernel(flags='fast-math')
    def _run_cleanup(self):
        """
        Set all devices back to their original state.
        """
        self.core.break_realtime()
        self.urukul0_cpld.set_profile(0)
        self.core.break_realtime()
        self.urukul1_cpld.set_profile(0)
        self.core.break_realtime()

    @rpc(flags='async')
    def update_dataset(self, *args):
        """
        Records data from the main sequence in the experiment dataset.

        Parameters passed to this function will be converted into a 1D array and added to the dataset.
        For efficiency, data is added by mutating indices of a preallocated dataset.
        Contains an internal iterator to keep track of the current index.
        """
        self.results[self._result_iter] = array(args)
        self._result_iter += 1


    # RUN - USER FUNCTIONS
    def run_initialize(self):
        """
        To be subclassed.

        todo: document
        """
        pass

    def loop_prepare(self):
        """
        To be subclassed.

        todo: document
        """
        pass

    def loop_run(self):
        """
        To be subclassed.

        todo: document
        """
        pass


    '''
    HasEnvironment Extensions
    '''

    def setattr_argument(self, *args, **kwargs):
        super().setattr_argument(*args, **kwargs)

        # add argument to _build_arguments (will grab after prepare)
        key, processor = args
        self._build_arguments[key] = None