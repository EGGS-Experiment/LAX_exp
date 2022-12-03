from artiq.experiment import *

import logging
from abc import ABC, abstractmethod
from inspect import getmembers, ismethod


logger = logging.getLogger("artiq.master.experiments")


class LAXExperiment(HasEnvironment, ABC):
    """
    Base class for LAX experiments.
    """

    kernel_invariants =     set()

    DDS_BOARD =             None
    DDS_CHANNEL =           None

    frequencies =           []
    amplitudes =            []


    # SETUP
    def build(self):
        """
        #todo document
        """
        # get core device
        self.setattr_device("core")

        # get dds channel
        urukul_cpld = self.get_dataset(self.DDS_BOARD, archive=False)
        urukul_chan = self.get_dataset(self.DDS_CHANNEL, archive=False)
        setattr(
            self,
            "dev",
            self.get_device('urukul{:d}_ch{:d}'.format(urukul_cpld, urukul_chan))
        )

        # get parameters from master dataset
        self._build_set_parameters()

        # set waveforms
        self._build_set_profiles()

        # steal all relevant methods of underlying TTL counter object so users
        # can directly call methods from this wrapper
        isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")
        device_functions = getmembers(self.dev, isDeviceFunction)
        for (function_name, function_object) in device_functions:
            setattr(self, function_name, function_object)


    def _build_set_parameters(self):
        """
        set parameters
        # todo document
        """
        # get frequencies
        for freq_param in self.frequencies:

            # reformat parameter name
            freq_param_new = freq_param.split('.')[-1]
            freq_param_new = freq_param_new.replace('mhz', 'ftw')

            # get parameter from dataset manager
            freq_val_mhz = self._HasEnvironment__dataset_mgr.ddb.get(freq_param)

            # set as parameter in dataset manager and HDF5 file
            self.setattr_dataset(
                freq_param_new,
                self.dev.frequency_to_ftw(freq_val_mhz * MHz),
                archive=True
                # todo: make dataset manager store parameters differently
            )

            # add parameter to kernel invariants
            self.kernel_invariants.add(freq_param_new)


        # get amplitudes
        for ampl_param in self.amplitudes:

            # reformat parameter name
            ampl_param_new = ampl_param.split('.')[-1]
            ampl_param_new = ampl_param_new.replace('pct', 'asf')

            # get parameter from dataset manager
            ampl_val_pct = self._HasEnvironment__dataset_mgr.ddb.get(ampl_param)

            # set as parameter in dataset manager and HDF5 file
            self.setattr_dataset(
                ampl_param_new,
                self.dev.amplitude_to_asf(ampl_val_pct / 100),
                archive=True
                # todo: make dataset manager store parameters differently
            )

            # add parameter to kernel invariants
            self.kernel_invariants.add(ampl_param_new)

    @kernel
    def _build_set_profiles(self):
        """
        #todo: document
        """
        pass
