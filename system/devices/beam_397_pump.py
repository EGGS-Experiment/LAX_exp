from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam397Pump(LAXDevice):
    """
    Wrapper for the 397nm pump beam.
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "pump"

    parameters = {
        'freq_cooling_ftw':                 ('beams.freq_mhz.freq_pump_cooling_mhz',    mhz_to_ftw),
        'freq_readout_ftw':                 ('beams.freq_mhz.freq_pump_readout_mhz',    mhz_to_ftw),
        'freq_rescue_ftw':                  ('beams.freq_mhz.freq_pump_rescue_mhz',     mhz_to_ftw),

        'ampl_cooling_asf':                 ('beams.ampl_pct.ampl_pump_cooling_pct',    pct_to_asf),
        'ampl_readout_asf':                 ('beams.ampl_pct.ampl_pump_readout_pct',    pct_to_asf),
        'ampl_rescue_asf':                  ('beams.ampl_pct.ampl_pump_rescue_pct',     pct_to_asf)
    }
    core_devices = {
        'beam': 'urukul1_ch1'
    }

    @kernel(flags={"fast-math"})
    def prepare_device(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=1)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def cooling(self):
        # set cooling profile
        with parallel:
            self.beam.cpld.set_profile(0)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def readout(self):
        # set readout profile
        with parallel:
            self.beam.cpld.set_profile(1)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def rescue(self):
        # set rescue profile
        with parallel:
            self.beam.cpld.set_profile(2)
            delay_mu(TIME_PROFILESWITCH_DELAY_MU)
