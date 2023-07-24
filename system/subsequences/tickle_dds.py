from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class TickleDDS(LAXSubsequence):
    """
    Subsequence: Tickle DDS

    Heat the ion by applying an RF signal from a DDS at the secular frequency.
    """
    name = 'tickle_dds'


    def build_subsequence(self):
        self.setattr_argument('time_tickle_us',         NumberValue(default=1000, ndecimals=3, step=10, min=1, max=1000000), group='tickle_dds')
        self.setattr_argument('att_tickle_db',          NumberValue(default=30, ndecimals=1, step=0.5, min=8, max=31.5), group='tickle_dds')

        # get relevant devices
        self.setattr_device('dds_modulation')

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.time_tickle_mu =                           self.core.seconds_to_mu(self.time_tickle_us * us)
        self.att_tickle_mu =                            att_to_mu(self.att_tickle_db * dB)

    @kernel(flags={"fast-math"})
    def run(self):
        # set dds attenuation here - ensures that dds channel will have correct attenuation
        self.dds_modulation.set_att_mu(self.att_tickle_mu)
        # change to standard profile
        self.dds_modulation.set_profile(0)
        # reset signal phase
        self.dds_modulation.reset_phase()

        # tickle for given time
        self.dds_modulation.on()
        delay_mu(self.time_tickle_mu)
        self.dds_modulation.off()
