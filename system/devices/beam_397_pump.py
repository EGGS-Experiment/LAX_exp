from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

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
        "sw", "cpld",

        "freq_cooling_ftw", "freq_readout_ftw", "freq_rescue_ftw",
        "ampl_cooling_asf", "ampl_readout_asf", "ampl_rescue_asf",
        "att_pump_mu",

        "profile_cooling", "profile_readout", "profile_rescue",
    }

    def prepare_device(self):
        # re-alias relevant base devices
        self.sw =   self.beam.sw
        self.cpld = self.beam.cpld

        # define profiles here for clarity/external use
        # todo: maybe move to build?
        self.profile_cooling =  0
        self.profile_readout =  1
        self.profile_rescue =   2

        # get attenuations
        # todo: check that attenuation is valid
        self.att_pump_mu =  self.get_parameter('att_pump_db', group='beams.att_db',
                                               override=False, conversion_function=att_to_mu)

        # get frequency parameters
        self.freq_cooling_ftw = self.get_parameter('freq_pump_cooling_mhz', group='beams.freq_mhz',
                                                   override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_readout_ftw = self.get_parameter('freq_pump_readout_mhz', group='beams.freq_mhz',
                                                   override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_rescue_ftw =  self.get_parameter('freq_pump_rescue_mhz', group='beams.freq_mhz',
                                                   override=False, conversion_function=hz_to_ftw, units=MHz)

        # get amplitude parameters
        self.ampl_cooling_asf = self.get_parameter('ampl_pump_cooling_pct', group='beams.ampl_pct',
                                                   override=False, conversion_function=pct_to_asf)
        self.ampl_readout_asf = self.get_parameter('ampl_pump_readout_pct', group='beams.ampl_pct',
                                                   override=False, conversion_function=pct_to_asf)
        self.ampl_rescue_asf =  self.get_parameter('ampl_pump_rescue_pct', group='beams.ampl_pct',
                                                   override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # get CPLD attenuations so we don't override them
        self.cpld.get_att_mu()
        self.core.break_realtime()

        # set waveforms for cooling, readout, and rescue
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=self.profile_cooling, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=self.profile_readout, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_rescue_ftw, asf=self.ampl_rescue_asf, profile=self.profile_rescue, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)

        # set attenuation
        self.set_att_mu(self.att_pump_mu)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        # set default profile on CPLD
        self.set_profile(DEFAULT_PROFILE)
        self.on()
        delay_mu(8000)

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        # enable RF switch onboard Urukul
        self.sw.on()
        delay_mu(8)

        # enable external RF switch
        self.rf_switch.off()
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        # disable RF switch onboard Urukul
        self.sw.off()
        delay_mu(8)

        # disable external RF switch
        self.rf_switch.on()
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def cooling(self) -> TNone:
        """
        Set doppler cooling profile.
        """
        self.cpld.set_profile(self.profile_cooling)
        self.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def readout(self) -> TNone:
        """
        Set readout (state-dependent fluorescence) profile.
        """
        self.cpld.set_profile(self.profile_readout)
        self.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def rescue(self) -> TNone:
        """
        Set rescue profile.
        """
        self.cpld.set_profile(self.profile_rescue)
        self.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        self.cpld.set_profile(profile_num)
        self.cpld.io_update.pulse_mu(8)

