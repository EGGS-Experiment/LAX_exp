from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


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

    @kernel(flags={"fast-math"})
    def prepare_device(self):
        # set base profile
        self.core.break_realtime()
        self.set_mu(mhz_to_ftw(1), asf=self.ampl_tickle_pct, profile=0)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)
