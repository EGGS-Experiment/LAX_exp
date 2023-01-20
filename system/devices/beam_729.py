from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam729(LAXDevice):
    """
    Device: Qubit Beam (729nm, polarized)

    Uses the DDS channel to drive the 729nm AOM in double-pass configuration.
    """
    name = "qubit"

    core_devices = {
        'beam': 'urukul0_ch1'
    }

    def prepare_device(self):
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False, conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self):
        # set waveform for carrier interrogation
        self.core.break_realtime()
        self.beam.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def on(self):
        self.beam.cfg_sw(True)

    @kernel(flags={"fast-math"})
    def off(self):
        self.beam.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def carrier(self):
        """
        Set carrier profile
        todo: document
        """
        self.beam.cpld.set_profile(0)
        self.beam.cpld.io_update.pulse_mu(8)
        delay_mu(TIME_PROFILESWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num):
        self.beam.cpld.set_profile(profile_num)
