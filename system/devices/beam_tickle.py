from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class BeamTickle(LAXDevice):
    """
    Device: DDS tickle beam

    Uses the DDS channel to apply a tickle on one of the radial
    """
    name = "tickle"

    core_devices = {
        'beam': 'urukul0_ch3'
    }

    def prepare_device(self):
        self.freq_tickle_ftw = hz_to_ftw(1 * MHz)
        self.ampl_tickle_pct = self.get_parameter('ampl_tickle_radial_pct', group='beams.ampl_pct', override=False, conversion=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        self.set_mu(self.freq_tickle_ftw, asf=self.ampl_tickle_pct, profile=0)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)
