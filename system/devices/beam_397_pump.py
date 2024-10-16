from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam397Pump(LAXDevice):
    """
    Device: Pump Beam (397nm)

    Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "pump"
    core_device = ('beam', 'urukul2_ch1')
    devices ={
        'rf_switch':    'ttl12'
    }
    kernel_invariants = {
        "freq_cooling_ftw",
        "freq_readout_ftw",
        "freq_rescue_ftw",
        "ampl_cooling_asf",
        "ampl_readout_asf",
        "ampl_rescue_asf"
    }

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
    def initialize_device(self) -> TNone:
        # set waveforms for cooling, readout, and rescue
        self.core.break_realtime()
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=1)
        self.core.break_realtime()
        self.set_mu(self.freq_rescue_ftw, asf=self.ampl_rescue_asf, profile=2)
        self.core.break_realtime()
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=3)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        with parallel:
            # enable RF switch onboard Urukul
            self.beam.sw.on()

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        with parallel:
            # disable RF switch onboard Urukul
            self.beam.sw.off()

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def cooling(self) -> TNone:
        """
        Set cooling profile
        todo: document
        """
        self.beam.cpld.set_profile(0)
        self.beam.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def readout(self) -> TNone:
        """
        Set readout profile
        todo: document
        """
        self.beam.cpld.set_profile(1)
        self.beam.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def rescue(self) -> TNone:
        """
        Set rescue profile
        todo: document
        """
        self.beam.cpld.set_profile(2)
        self.beam.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        self.beam.cpld.set_profile(profile_num)
        self.beam.cpld.io_update.pulse_mu(8)
