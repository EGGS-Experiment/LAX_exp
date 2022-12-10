from artiq.experiment import *
from inspect import getmembers, ismethod

from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam397Probe(LAXDevice):
    """
    Wrapper for the 397nm probe beam (polarized).
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "probe"

    parameters = {
        'freq_spinpol_ftw':          ('beams.freq_mhz.freq_probe_spinpol_mhz',      mhz_to_ftw),
        'ampl_spinpol_asf':          ('beams.ampl_pct.ampl_probe_spinpol_pct',      pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch0'
    }

    def prepare_device(self):
        # list of cfg_sw functions to break out
        sw_functions = ['on', 'off', 'pulse', 'pulse_mu']

        # verifies that a function is not magic
        isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (func_obj.__name__ is not "__init__")

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
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=1)
