from artiq.experiment import *

# todo: set base device name

class Beam_Urukul(HasEnvironment):
    """
    A generic "beam" object based off an Urukul DDS channel.
    """

    DDS_NAME = None
    frequencies = []
    amplitudes = []
    kernel_invariants = {}


    # BUILD
    def build(self):
        """
        #todo document
        """
        # get core devices
        self.setattr_device("core")
        self.dev = self.get_device(self.DDS_NAME)

        # get parameters from master dataset
        self._build_set_parameters()

        # set waveforms
        self._build_set_profiles()

    def _build_set_parameters(self):
        """
        set parameters
        # todo document
        """
        # get frequencies
        for freq_param in self.frequencies:

            # reformat parameter name
            freq_param_new = freq_param.split('.')[-1]
            freq_param_new.replace('mhz', 'ftw')

            # get parameter from dataset manager
            freq_val_mhz = self.__dataset_mgr.ddb.get(freq_param)

            # set as parameter in dataset manager and HDF5 file
            self.__dataset_mgr.set(
                freq_param_new,
                self.dev.frequency_to_ftw(freq_val_mhz),
                archive=False
                # todo: make dataset manager store parameters differently
            )

            # add parameter to kernel invariants
            self.kernel_invariants |= {freq_param}

        # get amplitudes
        for ampl_param in self.amplitudes:

            # reformat parameter name
            ampl_param_new = ampl_param.split('.')[-1]
            ampl_param_new.replace('pct', 'asf')

            # get parameter from dataset manager
            ampl_val_pct = self.__dataset_mgr.ddb.get(ampl_param)

            # set as parameter in dataset manager and HDF5 file
            self.__dataset_mgr.set(
                ampl_param_new,
                self.dev.amplitude_to_asf(ampl_val_pct),
                archive=False
                # todo: make dataset manager store parameters differently
            )

            # add parameter to kernel invariants
            self.kernel_invariants |= {ampl_param}

    @kernel
    def _build_set_profiles(self):
        """
        #todo: document
        """
        pass
