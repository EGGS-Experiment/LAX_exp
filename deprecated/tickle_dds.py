from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class TickleDDS(LAXSubsequence):
    """
    Subsequence: Tickle DDS

    Heat the ion by applying an RF signal from a DDS at the secular frequency.
    """
    name = 'tickle_dds'
    kernel_invariants = {
        "time_tickle_mu", "att_tickle_mu"
    }

    def build_subsequence(self):
        self.setattr_argument('time_tickle_us', NumberValue(default=100, precision=3, step=100, min=1, max=1000000),
                              group=self.name)
        self.setattr_argument('att_tickle_db',  NumberValue(default=10, precision=1, step=0.5, min=0, max=31.5),
                              group=self.name)

        # get relevant devices
        self.setattr_device('dds_dipole')

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.time_tickle_mu =   self.core.seconds_to_mu(self.time_tickle_us * us)
        self.att_tickle_mu =    att_to_mu(self.att_tickle_db * dB)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        # configure DDS here to ensure that correct attenuation during DMA recording
        self.dds_dipole.set_att_mu(self.att_tickle_mu)
        self.dds_dipole.set_profile(0)
        self.dds_dipole.set_phase_absolute()
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset DDS phase
        self.dds_dipole.reset_phase()

        # tickle for given time
        self.dds_dipole.on()
        delay_mu(self.time_tickle_mu)
        self.dds_dipole.off()
