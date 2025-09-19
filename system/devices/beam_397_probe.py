from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam397Probe(LAXDevice):
    """
    Device: Probe Beam (397nm, polarized)

    Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "probe"
    core_device = ('beam', 'urukul2_ch0')
    kernel_invariants = {
        "cpld", "sw",
        "freq_spinpol_ftw", "freq_rescue_ftw", "ampl_spinpol_asf", "ampl_rescue_asf",
        "att_probe_mu",
    }

    def prepare_device(self):
        # re-alias relevant base devices
        self.sw =   self.beam.sw
        self.cpld = self.beam.cpld

        # get attenuations
        # todo: check that attenuation is valid
        self.att_probe_mu = self.get_parameter('att_probe_db', group='beams.att_db',
                                               override=False, conversion_function=att_to_mu)

        # get frequency parameters
        self.freq_spinpol_ftw = self.get_parameter('freq_probe_spinpol_mhz', group='beams.freq_mhz',
                                                   override=False, conversion_function=hz_to_ftw, units=MHz)
        self.freq_rescue_ftw =  self.get_parameter('freq_probe_rescue_mhz', group='beams.freq_mhz',
                                                   override=False, conversion_function=hz_to_ftw, units=MHz)

        # get amplitude parameters
        self.ampl_spinpol_asf = self.get_parameter('ampl_probe_spinpol_pct', group='beams.ampl_pct',
                                                   override=False, conversion_function=pct_to_asf)
        self.ampl_rescue_asf =  self.get_parameter('ampl_probe_rescue_pct', group='beams.ampl_pct',
                                                   override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # get CPLD attenuations so we don't override them
        self.cpld.get_att_mu()
        self.core.break_realtime()

        # set cooling and readout profiles
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=1, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_rescue_ftw, asf=self.ampl_rescue_asf, profile=2, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(8000)

        # set attenuation
        self.set_att_mu(self.att_probe_mu)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        self.sw.off()
        delay_mu(5000)

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        self.sw.on()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        self.sw.off()
        delay_mu(8)
