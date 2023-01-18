from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam397Probe(LAXDevice):
    """
    Device: 397nm probe beam (polarized)

    Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "probe"

    core_devices = {
        'beam': 'urukul1_ch0'
    }

    def prepare_device(self):
        self.freq_spinpol_ftw = self.get_parameter('freq_spinpol_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_spinpol_asf = self.get_parameter('ampl_spinpol_asf', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # set cooling and readout profiles
        # todo: tune delays
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=1)
        self.core.break_realtime()
        self.set_mu(self.freq_spinpol_ftw, asf=self.ampl_spinpol_asf, profile=2)

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)
