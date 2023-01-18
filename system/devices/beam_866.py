from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam866(LAXDevice):
    """
    Device: 866nm cooling repump

    Uses the DDS channel to drive an AOM.
    """
    name = "repump_cooling"

    core_devices = {
        'beam': 'urukul1_ch2'
    }

    def prepare_device(self):
        self.freq_repump_cooling_ftw = self.get_parameter('freq_repump_cooling_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_repump_cooling_asf = self.get_parameter('ampl_repump_cooling_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()
        self.beam.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=2)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)
