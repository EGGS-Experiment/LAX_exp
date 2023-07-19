from artiq.experiment import *

import os
import time
import h5py
import socket
import logging

from sipyco import pyon
from datetime import datetime
from numpy import array, zeros
from abc import ABC, abstractmethod

logger = logging.getLogger("artiq.master.experiments")

from LAX_exp.base import LAXEnvironment, LAXDevice, LAXSequence, LAXSubsequence
from LAX_exp.base.manager_wrappers import _write_to_group


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
        self.setattr_device('ttl20')
        self.setattr_device('ttl21')
        self.setattr_device('ttl22')
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
        ***todo*** prepare_device is called before prepare_experiment to ensure that any values in prepare_device
        _prepare_experiment is called after prepare_experiment since subsequence instantiation may require
            values computed only in prepare_experiment.
        """
        # store arguments in dataset manager
        self._save_arguments()

        # todo: call device prepare method first

        # call user-defined prepare function
        self.prepare_experiment()

        # collate initialize functions to speed up initialization
        # note: devices should be initialized first
        _initialize_device_list =               [obj
                                                for obj in self.children
                                                if isinstance(obj, LAXDevice)]
        _initialize_subsequence_list =          [obj
                                                for obj in self.children
                                                if isinstance(obj, LAXSubsequence)]
        _initialize_sequence_list =             [obj
                                                for obj in self.children
                                                if isinstance(obj, LAXSequence)]


        # write code which initializes the relevant modules entirely on the kernel, without any RPCs
        _initialize_code =                      "self.core.reset()\n"

        # code to initialize devices
        for i, obj in enumerate(_initialize_device_list):
            if isinstance(obj, LAXDevice):
                setattr(self, '_LAXDevice_{}'.format(i), obj)
                _initialize_code += "self._LAXDevice_{}.initialize_device()\n".format(i)

        # code to initialize subsequences
        for i, obj in enumerate(_initialize_subsequence_list):
            if isinstance(obj, LAXSubsequence):
                setattr(self, '_LAXSubsequence_{}'.format(i), obj)
                _initialize_code += "self._LAXSubsequence_{}.initialize_subsequence()\n".format(i)

        # code to initialize sequences
        for i, obj in enumerate(_initialize_sequence_list):
            if isinstance(obj, LAXSequence):
                setattr(self, '_LAXSequence_{}'.format(i), obj)
                _initialize_code += "self._LAXSequence_{}.initialize_sequence()\n".format(i)
        _initialize_code += "self.core.break_realtime()"

        # create kernel from code string and set as _initialize_experiment
        initialize_func = kernel_from_string(["self"], _initialize_code)
        setattr(self, '_initialize_experiment', initialize_func)


        # call prepare methods of all child objects
        self.call_child_method('prepare')

        # create dataset for results
        self.set_dataset('results', zeros(self.results_shape))
        self.setattr_dataset('results')

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

        # start counting initialization time
        time_init_start = datetime.timestamp(datetime.now())

        # call initialize_* functions for all children in order of abstraction level (lowest first)
        self._initialize_experiment(self)

        # record initialization time
        time_init_stop = datetime.timestamp(datetime.now())
        print('\tInitialize Time: {:.2f}'.format(time_init_stop - time_init_start))

        # call user-defined initialize function
        # todo: see if we can move this into _initialize_experiment for speed
        self.initialize_experiment()

        # get DMA handles for subsequences recorded onto DMA
        # todo: move _load_dma onto initialize for speed
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
        pass

    @kernel(flags={"fast-math"})
    def _run_cleanup(self):
        """
        Set all devices back to their original state.
        """
        self.core.reset()

        # reset hardware to allow use
        with parallel:

            # reset qubit board
            with sequential:
                self.urukul0_cpld.set_profile(0)
                self.urukul0_cpld.io_update.pulse_mu(8)
                self.urukul0_cpld.cfg_switches(0b0000)

            # reset main board to rescue parameters
            with sequential:
                self.urukul1_cpld.set_profile(0)
                self.urukul1_cpld.io_update.pulse_mu(8)
                self.urukul1_cpld.cfg_switches(0b1110)

            # enable all RF switches
            self.ttl20.off()
            self.ttl21.off()
            self.ttl22.off()

        self.core.break_realtime()

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

    @rpc(flags={"async"})
    def update_results(self, *args):
        """
        Records data from the main sequence in the experiment dataset.

        Parameters passed to this function will be converted into a 1D array and added to the dataset.
        For efficiency, data is added by mutating indices of a preallocated dataset.
        Contains an internal iterator to keep track of the current index.
        """
        self.mutate_dataset('results', self._result_iter, array(args))
        # self.set_dataset('management.completion_pct', self._result_iter / len(self.results), broadcast=True)
        self._result_iter += 1

    @property
    @abstractmethod
    def results_shape(self):
        """
        todo: document
        """
        pass

    def write_results(self, exp_params):
        """
        Write arguments, datasets, and parameters in a well-structured format
        that uses the capabilities of HDF5.
        """
        # set variables
        rid = exp_params["rid"]
        start_time = time.localtime(exp_params["start_time"])
        expid = exp_params["expid"]
        # todo: try to get default save dir list
        # try:
        #     th0 = self.get_dataset('management.dataset_save_locations')
        #     print(th0)
        # except Exception as e:
        #     print(e)
        #     print('whoops')
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

                    # save expid separately to allow dashboard to run from hdf5 file
                    # note: expid has already been converted to hdf5-savable form in main artiq package
                    f["expid"] = expid

                    # store experiment parameters in a separate group as attributes
                    experiment_group = f.create_group("experiment")
                    for k, v in exp_params.items():
                        _write_to_group(experiment_group.attrs, k, v)

                    # get system parameters from labrad and save as attributes in a separate group
                    sys_params = self._save_labrad()
                    system_group = f.create_group("system")
                    for k, v in sys_params.items():
                        _write_to_group(system_group.attrs, k, v)

            # catch any errors
            except Exception as e:
                # todo: create error message
                print("Warning: unable to create and save file in LAX format: {}".format(e))

    def _save_labrad(self):
        """
        Extract system parameters from LabRAD and format for saving with the experiment results.
        """
        # todo: get wavemeter snapshot
        try:

            # import relevant modules
            import labrad

            # get default labrad connection values
            LABRADHOST =       os.environ['LABRADHOST']
            LABRADPORT =       os.environ['LABRADPORT']
            LABRADPASSWORD =   os.environ['LABRADPASSWORD']
            # create a synchronous connection to labrad
            # cxn = labrad.connect(LABRADHOST, port=LABRADPORT, name="{:s}_({:s})".format("ARTIQ_EXP", gethostname()), username="", password=LABRADPASSWORD)
            cxn = labrad.connect(name="{:s}_({:s})".format("ARTIQ_EXP", socket.gethostname()), username="", password=LABRADPASSWORD)
            # cxn_wm = labrad.connect(name="{:s}_({:s})".format("ARTIQ_EXP", socket.gethostname()), username="", password=LABRADPASSWORD)

            # connect to relevant servers
            rf = cxn.rf_server
            dc = cxn.dc_server
            ls = cxn.lakeshore336_server

            # get rf values
            rf.select_device()
            rf_freq_hz = rf.frequency()
            rf_ampl_dbm = rf.amplitude()

            # import dc config and get relevant parameters
            from EGGS_labrad.config.dc_config import dc_config
            active_channels = dc_config.channeldict

            # extract dc values of channels in use
            dc_vals = {}
            for channel_name, channel_params in active_channels.items():

                try:
                    # get channel number
                    channelnum = channel_params['num']

                    # store channel number
                    key_channelnum = "dc_channel_num_{:s}".format(channel_name)
                    dc_vals[key_channelnum] = channelnum

                    # store channel voltage
                    key_voltage = "dc_voltage_{:s}".format(channel_name)
                    voltage_v = dc.voltage(channelnum)
                    dc_vals[key_voltage] = voltage_v

                except Exception as e:
                    pass

            # get temperature values
            temp_vals_k = ls.read_temperature()
            temp_vals_k = array(list(temp_vals_k))

            # store all system values in a combined dict
            sys_vals = {
                "rf_freq_hz":   rf_freq_hz,
                "rf_ampl_dbm":  rf_ampl_dbm,
                "temp_vals_k":  temp_vals_k
            }
            sys_vals.update(dc_vals)

            # return combined system values
            return sys_vals


        except Exception as e:
            pass

        finally:
            # ensure labrad connection disconnects
            cxn.disconnect()

        # return empty dict to handle exception case
        return {}
