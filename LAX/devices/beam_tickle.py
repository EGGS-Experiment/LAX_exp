from artiq.experiment import *
from inspect import getmembers, ismethod

from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class BeamTickle(LAXDevice):
    """
    Wrapper for the tickle beam.
        Uses the DDS channel to apply a tickle on one of the radial
    """
    name = "tickle"

    parameters = {
        'ampl_tickle_pct':          ('beams.ampl_pct.ampl_tickle_radial_pct',       pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul0_ch3'
    }

    @kernel(flags='fast-math')
    def prepare_hardware(self):
        # set base profile
        self.core.break_realtime()
        self.set_mu(mhz_to_ftw(1), asf=self.ampl_spinpol_asf, profile=0)

    def prepare_class(self):
        # list of functions to break out
        sw_functions = ['on', 'off', 'pulse', 'pulse_mu']

        # verifies that a function is not magic
        isDeviceFunction = lambda func_obj: (callable(func_obj)) and (ismethod(func_obj)) and (
                    func_obj.__name__ is not "__init__")

        # steal all relevant methods of underlying device objects so users can directly call methods from this wrapper
        for (function_name, function_object) in getmembers(self.beam.cfg_sw, isDeviceFunction):
            if function_name in sw_functions:
                setattr(self, function_name, function_object)
