from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: on/off, set carrier via DUC, pulse shape


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF
    """
    name = "phaser_eggs"
    core_device = ('beam', 'phaser0')

    def prepare_device(self):
        # get frequency parameters
        self.freq_cooling_ftw =     self.get_parameter('freq_pump_cooling_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_readout_ftw =     self.get_parameter('freq_pump_readout_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_rescue_ftw =      self.get_parameter('freq_pump_rescue_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)

        # get amplitude parameters
        self.ampl_cooling_asf =     self.get_parameter('ampl_pump_cooling_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.ampl_readout_asf =     self.get_parameter('ampl_pump_readout_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.ampl_rescue_asf =      self.get_parameter('ampl_pump_rescue_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # set waveforms for cooling, readout, and rescue
        self.core.break_realtime()
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=1)
        self.core.break_realtime()
        self.set_mu(self.freq_rescue_ftw, asf=self.ampl_rescue_asf, profile=2)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def on(self):
        with parallel:
            # enable RF switch onboard Urukul
            self.beam.cfg_sw(True)

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_RFSWITCH_DELAY_MU)


    @kernel(flags={"fast-math"})
    def off(self):
        with parallel:
            # disable RF switch onboard Urukul
            self.beam.cfg_sw(False)

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_RFSWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def cooling(self):
        """
        Set cooling profile
        todo: document
        """
        self.beam.cpld.set_profile(0)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def readout(self):
        """
        Set readout profile
        todo: document
        """
        self.beam.cpld.set_profile(1)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def rescue(self):
        """
        Set rescue profile
        todo: document
        """
        self.beam.cpld.set_profile(2)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)