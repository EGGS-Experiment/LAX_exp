from artiq.experiment import *
from inspect import getmembers, ismethod

from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam854(LAXDevice):
    """
    Wrapper for the 854nm qubit repump.
        Uses the DDS channel to drive an AOM.
    """
    name = "qubit_repump"

    parameters = {
        'freq_repump_qubit_ftw':            ('beams.freq_mhz.freq_repump_qubit_mhz',        mhz_to_ftw),
        'ampl_repump_qubit_asf':            ('beams.ampl_pct.ampl_repump_qubit_pct',        pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch3'
    }

    def prepare_device(self):
        # list of cfg_sw functions to break out
        sw_functions = ['on', 'off', 'pulse', 'pulse_mu']

        # verifies that a function is not magic
        isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (
                    func_obj.__name__ is not "__init__")

        # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
        for (function_name, function_object) in getmembers(self.beam.cfg_sw, isDeviceFunction):
            if function_name in sw_functions:
                setattr(self, function_name, function_object)

        # prepare dds profiles
        self.prepare_hardware()

    @kernel(flags='fast-math')
    def prepare_hardware(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
