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
        self.setattr_device('dds_modulation')

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.time_tickle_mu =   self.get_parameter('time_tickle_ms', override=True, conversion_function=seconds_to_mu, units=ms)
        self.att_tickle_mu =    self.get_parameter('att_tickle_db', override=True, conversion_function=seconds_to_mu, units=ms)

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
        delay_mu(self.time_qlms_heating_mu)
        self.dds_modulation.off()
