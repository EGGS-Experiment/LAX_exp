from artiq.experiment import *

import os
import time
import h5py
import logging
from abc import ABC
from numpy import array

logger = logging.getLogger("artiq.master.experiments")

from LAX_exp.base import LAXEnvironment, LAXDevice, LAXSequence, LAXSubsequence

# tmp remove
from datetime import datetime


class LAXExperiment(LAXEnvironment, ABC):
    """
    Base class for experiment objects.

    Runs a sequence a given number of times and records corresponding data.
    Instantiates and initalizes all relevant devices, sequences, and subsequences.

    Attributes:
        name                        str                     : the name of the sequence (must be unique). Will also be used as the core_dma handle.
    """
    # Class attributes
    name =                  None
    _dma_count =            0


    '''
    BUILD
    '''

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

    def prepare(self):
        """
        General construction of the experiment object.

        Called right before run (unlike build, which is called upon instantiation).
        _prepare_experiment is called after prepare_experiment since subsequence instantiation may require
            values computed only in prepare_experiment.
        """
        # store arguments in dataset manager
        self._save_arguments()

        # call user-defined prepare function
        self.prepare_experiment()

        # collate initialize functions to speed up initialization
        # note: devices should be initialized first
        _initialize_device_list =               [inst_name
                                                for (inst_name, inst_obj) in self.__dict__.items()
                                                if isinstance(inst_obj, LAXDevice)]
        _initialize_subsequence_list =          [inst_name
                                                for (inst_name, inst_obj) in self.__dict__.items()
                                                if isinstance(inst_obj, LAXSubsequence)]
        _initialize_sequence_list =             [inst_name
                                                for (inst_name, inst_obj) in self.__dict__.items()
                                                if isinstance(inst_obj, LAXSequence)]


        # write code which initializes the relevant modules entirely on the kernel, without any RPCs
        _initialize_code =                      "self.core.reset()\n"
        for device_name in _initialize_device_list:
            self.setattr_device(device_name)
            _initialize_code +=                 "self.{}.initialize_device()\n".format(device_name)
            _initialize_code +=                 "self.core.break_realtime()\n"
        for subsequence_name in _initialize_subsequence_list:
            _initialize_code +=                 "self.{}.initialize_subsequence()\n".format(subsequence_name)
            _initialize_code +=                 "self.core.break_realtime()\n"
        for sequence_name in _initialize_sequence_list:
            _initialize_code +=                 "self.{}.initialize_sequence()\n".format(sequence_name)
            _initialize_code +=                 "self.core.break_realtime()\n"
        _initialize_code += "self.core.break_realtime()"

        # create kernel from code string and set as _initialize_experiment
        initialize_func = kernel_from_string(["self"], _initialize_code)
        setattr(self, '_initialize_experiment', initialize_func)

        # todo: get a labrad snapshot
        # need: trap rf amp/freq/locking, 6x dc voltages & on/off, temp, pressure
        # need: wavemeter frequencies, DDS attenuation, B-fields

        # todo: create dataset to hold results
        # todo: maybe specify dimensionality of results
        #self.set_dataset('results', list())
        #self.setattr_dataset('results')


        # call prepare methods of all child objects
        self.call_child_method('prepare')

    def prepare_experiment(self):
        """
        To be subclassed.

        Used to pre-compute data, modify arguments, and further instantiate objects.
        """
        pass


    '''
    RUN
    '''

    def run(self):
        """
        Main sequence of the experiment.
        Repeat a given sequence a number of times.
        """
        # set up completion monitor
        self.set_dataset('management.completion_pct', 0., broadcast=True, persist=True, archive=False)

        # tmp remove
        time1 = datetime.timestamp(datetime.now())

        # initialize children
        self._initialize_experiment(self)

        # tmp remove
        time2 = datetime.timestamp(datetime.now())
        print('\tinitialize time: {:.2f}'.format(time2-time1))

        # call user-defined initialize function
        self.initialize_experiment()

        # get DMA handles for subsequences recorded onto DMA
        self.call_child_method('_load_dma')

        # run the main part of the experiment
        self.run_main()

        # set devices back to their default state
        self._run_cleanup()

    @kernel(flags={"fast-math"})
    def _initialize_experiment(self):
        """
        Call the initialize functions of devices and sub/sequences (in that order).
        """
        # todo: see if necessary
        pass

    @kernel(flags={"fast-math"})
    def _run_cleanup(self):
        """
        Set all devices back to their original state.
        """
        self.core.reset()

        # reset qubit board
        self.urukul0_cpld.set_profile(0)
        self.urukul0_cpld.io_update.pulse_mu(8)
        self.urukul0_cpld.cfg_switches(0b0000)
        self.core.break_realtime()

        # reset main board to rescue parameters
        self.urukul1_cpld.set_profile(2)
        self.urukul1_cpld.io_update.pulse_mu(8)
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

    def initialize_experiment(self):
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
        # todo: repetitions sort out
        for trial_num in range(self.repetitions):

            # run the trial
            try:
                self.run_loop()

            # allow clean termination
            except TerminationRequested:
                break

    def run_loop(self):
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
                    self._LAXEnvironment__dataset_mgr.write_hdf5(f)

                    # store experiment parameters in a separate group as attributes
                    experiment_group = f.create_group("experiment")
                    for k, v in exp_params.items():
                        experiment_group.attrs[k] = v

            # catch any errors
            except Exception as e:
                print("Warning: unable to create and save file in LAX format: {}".format(e))
