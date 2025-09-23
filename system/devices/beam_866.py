from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam866(LAXDevice):
    """
    Device: Repump Cooling Beam (866nm)

    Uses the DDS channel to drive an AOM.
    """
    name = "repump_cooling"
    core_device = ('beam', 'urukul2_ch2')
    kernel_invariants = {
        "cpld", "sw",
        "freq_repump_cooling_ftw", "ampl_repump_cooling_asf", "att_repump_cooling_mu",
    }

    def prepare_device(self):
        # re-alias relevant base devices
        self.sw =   self.beam.sw
        self.cpld = self.beam.cpld

        # get attenuations
        # todo: check that attenuation is valid
        self.att_repump_cooling_mu =    self.get_parameter('att_repump_cooling_db', group='beams.att_db',
                                                           override=False, conversion_function=att_to_mu)

        # get beam parameters
        self.freq_repump_cooling_ftw = self.get_parameter('freq_repump_cooling_mhz', group='beams.freq_mhz',
                                                          override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_repump_cooling_asf = self.get_parameter('ampl_repump_cooling_pct', group='beams.ampl_pct',
                                                          override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # get CPLD attenuations so we don't override them
        self.cpld.get_att_mu()
        self.core.break_realtime()

        # set relevant beam parameters
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=2, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)

        # set attenuation
        self.set_att_mu(self.att_repump_cooling_mu)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        self.sw.on()
        delay_mu(5000)

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        self.sw.on()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        self.sw.off()
        delay_mu(8)
