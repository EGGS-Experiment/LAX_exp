from artiq.experiment import *

import os
import time
import h5py
import logging
from abc import ABC
from numpy import array

logger = logging.getLogger("artiq.master.experiments")

from LAX_exp.base import LAXEnvironment


class LAXExperiment(LAXEnvironment, ABC):
    """
    Base class for experiment objects.

    Runs a sequence a given number of times and records corresponding data.
    Instantiates and initalizes all relevant devices, sequences, and subsequences.

    Attributes:
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
    """
    # Class attributes
    _dma_count = 0

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
        setattr(self,'_result_iter', 0)

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
        self.save_arguments()
        self._prepare_experiment()
        self.prepare_experiment()
        self.call_child_method('prepare')

    def _prepare_experiment(self):
        """
        General construction of the experiment object.
        Must happen after the user-defined prepare_experiment method.
        """
        # todo: get a labrad snapshot
        # need: trap rf amp/freq/locking, 6x dc voltages & on/off, temp, pressure
        # need: wavemeter frequencies
        # need: DDS attenuation
        # need: B-fields

        # create dataset to hold results
        #self.set_dataset('results', list())
        #self.setattr_dataset('results')
        pass

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
        # todo: completion monitor stuff
        self.set_dataset('management.completion_pct', 0., broadcast=True, persist=True, archive=False)

        # set up the run by initializing devices, etc.
        self.run_initialize()

        # get DMA handles for subsequences recorded onto DMA
        self.call_child_method('_load_dma')

        # run the main part of the experiment
        self.run_main()

        # set devices back to their default state
        self._run_cleanup()

    @kernel(flags={"fast-math"})
    def _run_cleanup(self):
        """
        Set all devices back to their original state.
        """
        self.core.break_realtime()

        # reset qubit board
        self.urukul0_cpld.set_profile(0)
        self.urukul0_cpld.cfg_switches(0b0000)
        self.core.break_realtime()

        # reset main board to rescue
        self.urukul1_cpld.set_profile(2)
        self.urukul1_cpld.cfg_switches(0b1110)
        self.core.break_realtime()

    @rpc(flags={"async"})
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

    def run_main(self):
        """
        Can be subclassed.

        todo: document
        """
        # repeat the experiment a given number of times
        for trial_num in range(self.repetitions):

            # run the trial
            try:
                self.loop()

            # allow clean termination
            except TerminationRequested:
                break

    def loop(self):
        """
        To be subclassed.

        todo: document
        """
        pass


    '''
    Results
    '''

    def write_results(self, exp_params):
        """
        Write arguments, datasets, and parameters in a well-structured format
        that uses the capabilities of HDF5.
        """
        # set variables
        rid = exp_params["rid"]
        start_time = time.localtime(exp_params["start_time"])
        save_dir_list = ['Z:\\Motion\\Data']

        # save to all relevant directories
        for save_dir in save_dir_list:

            try:
                # format file name and save directory
                filedir = os.path.join(
                    save_dir,
                    time.strftime("%Y-%m-%d", start_time)
                )
                filename = "{:09}-{}.h5".format(rid, self.name)

                # enter save directory
                os.makedirs(filedir, exist_ok=True)
                os.chdir(filedir)

                # write data
                with h5py.File(filename, "w") as f:

                    # save data from experiment via the dataset manager of the LAXExperiment
                    self.__dataset_mgr.write_hdf5(f)

                    # store experiment parameters in a separate group as attributes
                    experiment_group = f.create_group("experiment")
                    for k, v in exp_params.items():
                        experiment_group.attrs[k] = v

            # catch any errors
            except Exception as e:
                print("Warning: unable to create and save file in LAX format: {}".format(e))
