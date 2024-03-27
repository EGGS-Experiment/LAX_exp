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
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # management
        self.setattr_device("scheduler")
        self.setattr_device("ccb")

        # science-related hardware devices
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_cpld')
        self.setattr_device('urukul2_cpld')

        self.setattr_device('urukul0_ch0')
        self.setattr_device('urukul1_ch1')
        self.setattr_device('urukul1_ch3')
        self.setattr_device('urukul2_ch2')
        self.setattr_device('urukul2_ch3')

        self.setattr_device('ttl12')
        self.setattr_device('ttl13')
        self.setattr_device('ttl14')

        self.setattr_device('phaser0')

        # set looping iterator for the _update_results method
        setattr(self, '_result_iter', 0)

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

        # note: we call prepare_device and prepare_subsequence here
        # to ensure that any attributes instantiated here will be accessible by prepare_experiment
        self.call_child_method("prepare_device")
        self.call_child_method("prepare_subsequence")

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


        # todo: somehow call prepare method for everything not devices
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

        try:
            # run the main part of the experiment
            self.run_main()
        except TerminationRequested:
            # run cleanup, bypassing any analysis methods
            self._run_cleanup()
            print('\tExperiment successfully terminated.')
            raise TerminationRequested

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
        self.core.break_realtime()
        # todo: set attenuations for parametric and modulation DDSs

        # reset hardware to allow use by users
        with parallel:
            # reset qubit board
            with sequential:
                self.urukul0_cpld.set_profile(0)
                self.urukul0_cpld.cfg_switches(0b0000)

            # reset motional board to rescue parameters
            with sequential:
                self.urukul1_cpld.set_profile(0)
                self.urukul1_cpld.cfg_switches(0b0000)

            # reset main board to rescue parameters
            with sequential:
                self.urukul2_cpld.set_profile(0)
                self.urukul2_cpld.cfg_switches(0b1110)

        delay_mu(100)
        with sequential:
            # enable all external RF switches
            with parallel:
                self.ttl12.off()
                self.ttl13.off()
                self.ttl14.off()

            # set urukul ttl switches to allow front-end access
            # since output is logical OR of the TTL state as well as the urukul configuration register state
            delay_mu(10)
            with parallel:
                self.urukul0_ch0.sw.off()  # 729nm
                self.urukul1_ch1.sw.off()  # parametric
                self.urukul1_ch3.sw.off()  # dipole
                self.urukul2_ch2.sw.off()  # 866nm
                self.urukul2_ch3.sw.off()  # 854nm

        # reset phaser attenuators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(40)
        self.phaser0.channel[1].set_att(31.5 * dB)

        # reset phaser oscillators
        for i in range(5):
            # synchronize to frame
            at_mu(self.phaser0.get_next_frame_mu())

            # clear oscillator frequencies
            with parallel:
                self.phaser0.channel[0].oscillator[i].set_frequency(0.)
                self.phaser0.channel[1].oscillator[i].set_frequency(0.)
                delay_mu(40)

            # clear oscillator amplitudes
            with parallel:
                self.phaser0.channel[0].oscillator[i].set_amplitude_phase(amplitude=0.)
                self.phaser0.channel[1].oscillator[i].set_amplitude_phase(amplitude=0.)
                delay_mu(40)

            # add slack
            self.core.break_realtime()

        # ensure all events in the FIFOs are completed before
        # we exit the kernel
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())

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

    @rpc
    def check_termination(self):
        """
        Check whether termination of the experiment has been requested.
        """
        if self.scheduler.check_termination():
            raise TerminationRequested

    @kernel(flags={"fast-math"})
    def noop(self):
        """
        A hardware no-op function to allow for customizable pulse sequence configuration.
        """
        delay_mu(0)


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

        # add error handling for experiment analysis
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

    def analyze_experiment(self):
        """
        To be subclassed.

        Used to process/analyze experiment results.
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
        self.set_dataset('management.completion_pct', round(100. * self._result_iter / len(self.results), 3), broadcast=True, persist=True, archive=False)
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

        # todo: try to get default save dir list
        # try:
        #     th0 = self.get_dataset('management.dataset_save_locations')
        #     print(th0)
        # except Exception as e:
        #     print(e)
        #     print('whoops')
        save_dir_list = [
            'Z:\\Motion\\Data'
        ]

        # save to all relevant directories
        for save_dir in save_dir_list:

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
            print("Warning: unable to retrieve and store trap RF values in dataset.")

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
