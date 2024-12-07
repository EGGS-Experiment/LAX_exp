from artiq.experiment import *

import os
import time
import h5py
import socket
import logging
import traceback

from datetime import datetime
from numpy import array, zeros, int32, nan
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
    name =          None
    _dma_count =    0


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
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # management
        self.setattr_device("scheduler")
        self.setattr_device("ccb")

        # set looping iterators for the _update_results method
        setattr(self, '_result_iter', 0)
        setattr(self, '_counts_iter', 0)

        # labrad termination checking
        self._termination_status_labrad = False

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
        # todo: consider synchronization - is it necessary? where should it be?
        # todo: maybe we should move the compilation part to its own function?
        """
        '''BEGIN PREPARE'''
        # get kernel invariants (if any)
        kernel_invariants = getattr(self, "kernel_invariants", set())
        self.kernel_invariants = kernel_invariants

        # store arguments in dataset manager
        self._save_arguments()

        # note: we call prepare_device and prepare_subsequence here
        # to ensure that any attributes instantiated here will be accessible by prepare_experiment
        self.call_child_method("prepare_device")
        self.call_child_method("prepare_subsequence")
        self.call_child_method("prepare_sequence")

        # call user-defined prepare function
        self.prepare_experiment()

        # todo: somehow call prepare method for everything not devices
        # todo: maybe - we could make LAX classes NOT put prepare_<class> under prepare?
        #       there's no reason they need to have a "prepare" method since they're not an EnvExp
        #       and they're not being run by themselves
        # call prepare methods of all child objects
        # note: this will call the prepare_<class> methods AGAIN, i.e. any preparation
        # that happens in prepare_experiment will get reset
        self.call_child_method('prepare')

        # create data structures for results
        self.set_dataset('results', zeros(self.results_shape))
        self.setattr_dataset('results')

        # create data structures for dynamic updates
        self._dynamic_reduction_factor = self.get_dataset('management.dynamic_plot_reduction_factor',
                                                          default=10, archive=False)
        self.kernel_invariants.add("_dynamic_reduction_factor")
        # preprocess values for completion monitoring
        self._completion_iter_to_pct = 100. / len(self.results)
        self.kernel_invariants.add("_completion_iter_to_pct")

        # get list of directories that result files should be saved to
        save_dir_list = {
            'Z:\\Motion\\Data',     # save to motion drive
            'D:\\Results'           # save to local backup drive
        }
        try:
            save_dir_list = save_dir_list | set(self.get_dataset('management.dataset_save_locations'))
        except Exception as e:
            print(repr(e))
        # set save_dir_list as an instance attribute
        self.save_dir_list = save_dir_list
        self.kernel_invariants.add("save_dir_list")


        '''SUBSCRIBE TO LABRAD FOR WARNINGS'''
        # get monitoring status from dataset manager
        monitor_status = self.get_dataset('management.safe_mode', default=False)
        if monitor_status:
            self.labrad_subscribe()


        '''COMPILE INITIALIZATION SEQUENCE'''
        # collate LAX objects to speed up experiment compilation and execution
        # note: objects should be initialized in the following order: devices, subsequences, sequences, then experiment
        _compile_device_list =              [obj
                                            for obj in self.children
                                            if isinstance(obj, LAXDevice)]
        _compile_subsequence_list =         [obj
                                            for obj in self.children
                                            if isinstance(obj, LAXSubsequence)]
        _compile_sequence_list =            [obj
                                            for obj in self.children
                                            if isinstance(obj, LAXSequence)]

        # write code which initializes the relevant modules entirely on the kernel, without any RPCs
        _initialize_code =          "self.core.reset()\n"
        # write code to get DMA handles for LAXSubsequences
        _initialize_load_DMA_code = "self.core.break_realtime()\n"

        # todo: should we set e.g. "_LAX<Object>_{}" as a kernel_invariant? or otherwise get their names?
        # todo: convert all of these for loops into list comprehensions
        # create code to initialize devices
        for i, obj in enumerate(_compile_device_list):
            if isinstance(obj, LAXDevice):
                setattr(self, '_LAXDevice_{}'.format(i), obj)
                self.kernel_invariants.add("_LAXDevice_{}".format(i))
                _initialize_code += "self._LAXDevice_{}.initialize_device()\n".format(i)

        # create code to initialize subsequences
        for i, obj in enumerate(_compile_subsequence_list):
            if isinstance(obj, LAXSubsequence):
                setattr(self, '_LAXSubsequence_{}'.format(i), obj)
                self.kernel_invariants.add("_LAXSubsequence_{}".format(i))
                _initialize_code += "self._LAXSubsequence_{}.initialize_subsequence()\n".format(i)
                # for LAXSubsequences only: get DMA handle for sequences recorded onto DMA
                _initialize_load_DMA_code += "self._LAXSubsequence_{}._load_dma()\n".format(i)

        # create code to initialize sequences
        for i, obj in enumerate(_compile_sequence_list):
            if isinstance(obj, LAXSequence):
                setattr(self, '_LAXSequence_{}'.format(i), obj)
                self.kernel_invariants.add("_LAXSequence_{}".format(i))
                _initialize_code += "self._LAXSequence_{}.initialize_sequence()\n".format(i)

        # call user-defined initialize function for the experiment
        _initialize_code += "self.initialize_experiment()\n"
        _initialize_code += "self.core.break_realtime()\n"
        # record DMA sequences and get DMA handles to play subsequences
        _initialize_code += _initialize_load_DMA_code
        _initialize_code += "self.core.break_realtime()"

        # create kernel from code string and set as _initialize_experiment
        initialize_func = kernel_from_string(["self"], _initialize_code)
        setattr(self, '_initialize_experiment', initialize_func)


        '''COMPILE CLEANUP SEQUENCE'''
        # collate cleanup functions to speed up experiment compilation and execution
        # note: objects should be cleaned up in the following order: experiments, sequences, subsequences, then devices
        # i.e. reverse order of initialization
        # write code which cleans up the relevant modules
        _cleanup_code = "self.core.break_realtime()\n"

        # call user-defined cleanup function for the experiment
        _cleanup_code += "self.cleanup_experiment()\n"
        _cleanup_code += "self.core.break_realtime()\n"

        # code to cleanup sequences
        for i, obj in enumerate(_compile_sequence_list):
            if isinstance(obj, LAXSequence):
                _cleanup_code += "self._LAXSequence_{}.cleanup_sequence()\n".format(i)

        # code to cleanup subsequences
        for i, obj in enumerate(_compile_subsequence_list):
            if isinstance(obj, LAXSubsequence):
                _cleanup_code += "self._LAXSubsequence_{}.cleanup_subsequence()\n".format(i)

        # code to cleanup devices
        for i, obj in enumerate(_compile_device_list):
            if isinstance(obj, LAXDevice):
                _cleanup_code += "self._LAXDevice_{}.cleanup_device()\n".format(i)

        # create kernel from code string and set as _cleanup_experiment
        _cleanup_code += "self.core.break_realtime()"
        cleanup_func = kernel_from_string(["self"], _cleanup_code)
        setattr(self, '_cleanup_experiment', cleanup_func)

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
        # set up dynamic datasets
        # note: this has to happen during run, otherwise we will overwrite other count plotters
        self.set_dataset('management.dynamic.completion_pct', 0., broadcast=True, persist=True, archive=False)
        # downsample counts for dynamic plotting
        dynamic_counts_len = (self.results_shape[0] // self._dynamic_reduction_factor) + 1
        dynamic_counts_arr = zeros(dynamic_counts_len, dtype=int32) * nan
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        dynamic_counts_arr[0] = 0
        self.set_dataset('temp.counts.trace', dynamic_counts_arr,
                         broadcast=True, persist=False, archive=False)

        # start counting initialization time
        time_global_start = datetime.timestamp(datetime.now())

        # call initialize_* functions for all children in order of abstraction level (lowest first)
        # note: this also calls the experiment's "initialize_experiment" function
        # note: needs to be passed "self" as an argument since it's a kernel_from_string
        self._initialize_experiment(self)

        # record initialization time
        time_init_stop = datetime.timestamp(datetime.now())
        time_expinit = time_init_stop - time_global_start
        self.set_dataset('time_init', time_expinit)
        print('\tInitialize Time:\t{:.2f}'.format(time_expinit))

        # run the main part of the experiment
        try:
            self.run_main()
        except TerminationRequested:
            print('\tExperiment successfully terminated.')
        except Exception as e:
            print('\tError during experiment: {}'.format(e))
            print(traceback.format_exc())
        finally:
            pass

        # clean up system/hardware to ensure system is left in a safe state
        # note: needs to be passed "self" as an argument since it's a kernel_from_string
        self._cleanup_experiment(self)

        # record total runtime
        time_run_stop = datetime.timestamp(datetime.now())
        time_exprun = time_run_stop - time_init_stop
        print('\tRun Time:\t\t{:.2f}'.format(time_exprun))
        self.set_dataset('time_run', time_exprun)

    @kernel(flags={"fast-math"})
    def _initialize_experiment(self) -> TNone:
        """
        Call the initialize functions of devices and sub/sequences (in that order).
        """
        pass

    def initialize_experiment(self) -> TNone:
        """
        To be subclassed.

        todo: document
        """
        pass

    def run_main(self) -> TNone:
        """
        Can be subclassed.

        todo: document
        """
        # repeat the experiment a given number of times
        # todo: repetitions sort out - who sets the argument for the exp?
        for trial_num in range(self.repetitions):
            # run the trial
            self.run_loop()

            # check if user requests experiment termination
            self.check_termination()

    def run_loop(self):
        """
        To be subclassed.

        todo: document
        """
        pass

    @kernel(flags={"fast-math"})
    def noop(self) -> TNone:
        """
        A hardware no-op function to allow for customizable pulse sequence configuration.
        Provided for convenience.
        """
        delay_mu(0)

    @kernel(flags={"fast-math"})
    def _cleanup_experiment(self) -> TNone:
        """
        Call the cleanup functions of sequences, subsequence, and devices and (in that order).
        """
        pass

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        """
        To be subclassed.
        todo: document
        """
        pass


    '''
    STATUS MONITORING
    '''

    @rpc
    def check_termination(self) -> TNone:
        """
        Check whether termination of the experiment has been requested.
        """
        if self.scheduler.check_termination() or self._termination_status_labrad:
            if self._termination_status_labrad:
                print("Critical experiment failure. Stopping experiment & cancelling all experiments.")
                self.cancel_all_experiments()
            raise TerminationRequested

    @rpc
    def labrad_subscribe(self) -> TNone:
        """
        Subscribe to local labrad warning server to get notifications
        when critical errors occur.
        """
        try:
            # import here to prevent repository scan issues
            from labrad.thread import startReactor
            from labrad.wrappers import connectAsync
            from twisted.internet.defer import inlineCallbacks, Deferred

            # start labrad's twisted reactor
            startReactor()
            d = Deferred()

            @inlineCallbacks
            def create_connection(msg):
                self.cxn_async = yield connectAsync(
                    os.environ["LABRADHOST"], port=7682, name="{:s} ({:s})".format("ARTIQ_EXP", socket.gethostname()),
                    username="", password=os.environ["LABRADPASSWORD"]
                )
                self.ws_async = self.cxn_async.warning_server
                self.ws_async.signal__wavemeter_unlock(959781)
                self.ws_async.addListener(listener=self._update_labrad_warnings, source=None, ID=959781)

            # fire deferred to create connection
            d.addCallback(create_connection)
            d.callback("\tDEFERRED: FIRED")
        except Exception as e:
            print("Unable to connect to labrad for warnings: {}".format(e))

    def _update_labrad_warnings(self, c, signal) -> TNone:
        """
        Automatically process warning updates from LabRAD,
        :param c: labrad context
        :param signal: the warning message/data.
        """
        # print("\n\t\tWARNING - CH{:d} UNLOCKED: {:s}".format(*signal))
        self.set_dataset('management.dynamic.ion_status',
                         "ERROR: {:s} UNLOCKED".format(signal[1]),
                         broadcast=True)
        self._termination_status_labrad = True

    @rpc
    def cancel_all_experiments(self):
        """
        Terminate all experiments in our pipeline.
        To be used in emergency situations (e.g. wavemeter unlocked) where we need to halt all activity.
        """
        # get scheduler itinerary
        sched = self.scheduler.get_status()

        # get all experiments in our pipeline
        rid_list = [
            rid
            for rid, exp_dict in sched.items()
            if (rid != self.scheduler.rid) and (exp_dict['pipeline'] == self.scheduler.pipeline_name)
               and (exp_dict['status'] != "running")
        ]
        rid_list.reverse()

        # delete remaining experiments
        for rid in rid_list:
            self.scheduler.delete(rid)


    '''
    ANALYZE
    '''

    def analyze(self):
        """
        General experimental analysis.

        Used to process/analyze experiment results.
        We separate the analyze stage for experiments into "analyze_experiment."

        This allows all other lower-class modules to straighforwardly define their own "analyze" methods
        and have them called by the parent experiment, as well as for those modules to run their analysis
        in isolation.
        """
        self.call_child_method('analyze')

        # only run analyze_experiment if exp ran to completion
        if self._result_iter == len(self.results):
            try:
                # get return from analyze_experiment method
                res_processed = self.analyze_experiment()

                # if we get a valid return, assume it is the processed result
                # of the experiment, and try to save it in the hdf5 file
                if res_processed is not None:
                    self.set_dataset('results_processed', res_processed)
            except Exception as e:
                print('\tWarning: encountered error during analyze_experiment.')
                print('\t\t:{:}'.format(repr(e)))
        else:
            print("Experiment results dataset not filled. Skipping analyze_experiment.")

    def analyze_experiment(self):
        """
        To be subclassed.

        Used to process/analyze experiment results.
        # todo: document what happens to any returns
        """
        return None


    '''
    RESULTS & DATASETS
    '''

    @rpc(flags={"async"})
    def update_results(self, *args):
        """
        Records data from the main sequence in the experiment dataset.

        Parameters passed to this function will be converted into a 1D array and added to the dataset.
        For efficiency, data is added by mutating indices of a preallocated dataset.
        Contains an internal iterator to keep track of the current index.
        """
        # store results in main dataset
        self.mutate_dataset('results', self._result_iter, array(args))

        # do intermediate processing
        if (self._result_iter % self._dynamic_reduction_factor) == 0:
            # plot counts in real-time to monitor ion death
            self.mutate_dataset('temp.counts.trace', self._counts_iter, args[1])
            self._counts_iter += 1

            # monitor completion status
            self.set_dataset('management.dynamic.completion_pct', round(self._result_iter * self._completion_iter_to_pct, 3),
                             broadcast=True, persist=True, archive=False)

        # increment result iterator
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
        rid =           exp_params["rid"]
        start_time =    time.localtime(exp_params["start_time"])
        expid =         exp_params["expid"]

        # save to all relevant directories - these are retrieved & stored in "prepare"
        for save_dir in self.save_dir_list:

            try:
                # format file name and save directory
                filedir = os.path.join(
                    save_dir,
                    time.strftime("%Y-%m", start_time),
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
                print("Warning: unable to create and save file in LAX format: {}".format(repr(e)))

    def _save_labrad(self):
        """
        Extract system parameters from LabRAD and format for saving with the experiment results.
        """
        # create holding dict for system parameters
        sys_vals = {}

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

            # create a synchronous connection to wavemeter labrad
            from EGGS_labrad.config.multiplexerclient_config import multiplexer_config
            cxn_wm = labrad.connect(multiplexer_config.ip, name="{:s}_({:s})".format("ARTIQ_EXP", socket.gethostname()), username="", password=LABRADPASSWORD)

            # get relevant system values
            sys_vals.update(self._save_labrad_rf(cxn))
            sys_vals.update(self._save_labrad_dc(cxn))
            sys_vals.update(self._save_labrad_temp(cxn))
            sys_vals.update(self._save_labrad_wm(cxn_wm))
            sys_vals.update(self._save_labrad_bfield(cxn))

        except Exception as e:
            print("Warning: error retrieving and saving LabRAD values in dataset.")

        finally:
            # ensure labrad connections disconnect
            cxn.disconnect()
            cxn_wm.disconnect()

        # return collated relevant system values
        return sys_vals

    def _save_labrad_rf(self, cxn):
        """
        Extract system RF parameters from LabRAD.
        Arguments:
            cxn     labrad_cxn  : a labrad connection object.
        Returns:
                    dict        : a dict of relevant system values.
        """
        # create holding dict
        sys_vals_rf = {}
        try:
            # set up RF server cxn
            rf = cxn.rf_server
            rf.select_device()

            # get RF values
            rf_freq_hz = rf.frequency()
            rf_ampl_dbm = rf.amplitude()

            # store rf values in holding dict
            sys_vals_rf = {
                "rf_freq_hz": rf_freq_hz,
                "rf_ampl_dbm": rf_ampl_dbm
            }
        except Exception as e:
            pass
            # print("Warning: unable to retrieve and store trap RF values in dataset.")

        return sys_vals_rf

    def _save_labrad_dc(self, cxn):
        """
        Extract system DC parameters from LabRAD.
        Arguments:
            cxn     labrad_cxn  : a labrad connection object.
        Returns:
                    dict        : a dict of relevant system values.
        """
        # create holding dict
        sys_vals_dc = {}
        try:
            # import DC config and get relevant parameters
            from EGGS_labrad.config.dc_config import dc_config
            dc_channels = dc_config.channeldict

            # extract values for DC channels in use
            dc = cxn.dc_server
            for channel_name, channel_params in dc_channels.items():
                try:
                    # get channel number
                    channelnum = channel_params['num']

                    # store channel number
                    key_channelnum = "dc_channel_num_{:s}".format(channel_name)
                    sys_vals_dc[key_channelnum] = channelnum

                    # store channel voltage
                    key_voltage = "dc_voltage_{:s}_v".format(channel_name)
                    voltage_v = dc.voltage(channelnum)
                    sys_vals_dc[key_voltage] = voltage_v
                except Exception as e:
                    pass

        except Exception as e:
            print("Warning: unable to retrieve and store trap DC values in dataset.")

        return sys_vals_dc

    def _save_labrad_temp(self, cxn):
        """
        Extract system temperature/cryo parameters from LabRAD.
        Arguments:
            cxn     labrad_cxn  : a labrad connection object.

        Returns:
                    dict        : a dict of relevant system values.
        """
        # create holding dict
        sys_vals_temp = {}
        try:
            # get temperature values from lakeshore
            ls = cxn.lakeshore336_server
            temp_vals_k = ls.read_temperature()

            # store temperature values in holding dict
            sys_vals_temp = {"temp_vals_k":  array(list(temp_vals_k))}
        except Exception as e:
            print("Warning: unable to retrieve and store temperature values in dataset.")

        return sys_vals_temp

    def _save_labrad_wm(self, cxn):
        """
        Extract system wavemeter parameters from LabRAD.
        Arguments:
            cxn     labrad_cxn  : a labrad connection object.
        Returns:
                    dict        : a dict of relevant system values.
        """
        # create holding dict

        # extract frequencies of wavemeter channels in use
        sys_vals_wm = {}
        try:
            # import wavemeter config and get relevant parameters
            from EGGS_labrad.config.multiplexerclient_config import multiplexer_config
            wm_channels = multiplexer_config.channels

            # extract values for wavemeter channels in use
            wm = cxn.multiplexerserver
            for channel_name, channel_params in wm_channels.items():
                try:
                    # get channel number
                    channelnum = channel_params[0]

                    # store channel frequency
                    key_frequency = "wm_frequency_{:s}_thz".format(channel_name)
                    freq_thz = wm.get_frequency(channelnum)
                    sys_vals_wm[key_frequency] = freq_thz
                except Exception as e:
                    pass

        except Exception as e:
            print("Warning: unable to retrieve and store laser wavemeter values in dataset.")

        return sys_vals_wm

    def _save_labrad_bfield(self, cxn):
        """
        Extract system B-field parameters from LabRAD.
        Arguments:
            cxn     labrad_cxn  : a labrad connection object.
        Returns:
                    dict        : a dict of relevant system values.
        """
        # create holding dict
        sys_vals_bfield = {}
        try:
            # set up respective server cxns
            ke = cxn.keithley_2231a_server
            gpp = cxn.gpp3060_server

            # get B-field power supply current values
            # bfield_ke_volts = ke.measure_voltage(1)
            # bfield_ke_amps = ke.measure_current(1)
            bfield_gpp_volts = gpp.measure_voltage(2)
            bfield_gpp_amps = gpp.measure_current(2)

            # store rf values in holding dict
            sys_vals_bfield = {
                # "bfield_keithley_volts":  bfield_ke_volts,
                # "bfield_keithley_amps":  bfield_ke_amps,
                "bfield_gpp3060_volts": bfield_gpp_volts,
                "bfield_gpp3060_amps":  bfield_gpp_amps
            }

            # tmp remove
            # todo: either move into a different section, or rename "_save_..." methods to "_check_..."
            # verify that main b-field current is valid
            if abs(bfield_gpp_amps - 0.7) > 0.2:
                print("Warning: Main B-field current coils are not their typical value.")
                print("\tCurrent Value (Main coils): {:.3f}A".format(bfield_gpp_amps))
            # tmp remove
        except Exception as e:
            print("Warning: unable to retrieve and store B-field current values in dataset.")

        return sys_vals_bfield
