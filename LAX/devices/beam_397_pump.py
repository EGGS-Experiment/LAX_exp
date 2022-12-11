from artiq.experiment import *
from inspect import getmembers, ismethod

from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf, us_to_mu


class Beam397Pump(LAXDevice):
    """
    Wrapper for the 397nm pump beam.
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "pump"

    parameters = {
        'freq_cooling_ftw':                 ('beams.freq_mhz.freq_pump_cooling_mhz',    mhz_to_ftw),
        'ampl_cooling_asf':                 ('beams.ampl_pct.ampl_pump_cooling_pct',    pct_to_asf),
        'freq_readout_ftw':                 ('beams.freq_mhz.freq_pump_readout_mhz',    mhz_to_ftw),
        'ampl_readout_asf':                 ('beams.ampl_pct.ampl_pump_readout_pct',    pct_to_asf),
        'time_profileswitch_delay_mu':      ('timing.time_profileswitch_delay_us',      us_to_mu)
    }
    core_devices = {
        'beam': 'urukul1_ch1'
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
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=1)

    @kernel(flags='fast-math')
    def cooling(self):
        # set cooling profile
        with parallel:
            self.beam.cpld.set_profile(0)
            delay_mu(self.time_profileswitch_delay_mu)

    @kernel(flags='fast-math')
    def readout(self):
        # set readout profile
        with parallel:
            self.beam.cpld.set_profile(1)
            delay_mu(self.time_profileswitch_delay_mu)
