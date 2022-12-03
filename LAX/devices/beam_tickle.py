from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class BeamTickle(LAXDevice):
    """
    Wrapper for the tickle beam.
        Uses the DDS channel to apply a tickle on one of the radial
    """
    name = "tickle"

    device_names = {'beam': 'urukul0_ch3'}
    device_parameters = {
        'ampl_tickle_pct': ('beams.ampl_pct.ampl_probe_redist_pct', pct_to_asf)
    }

    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set base profile
        self.core.break_realtime()
        self.set_mu(mhz_to_ftw(1), asf=self.ampl_redist_asf, profile=0)
