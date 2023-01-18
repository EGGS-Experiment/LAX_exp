from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SpinPolarization(LAXSubsequence):
    """
    Subsequence: Spin Polarization

    Place the ion in the S-1/2 mj=-1/2 state using the polarized 397 probe beam.
    """
    name = 'spin_polarization'

    def build_subsequence(self):
        self.setattr_device('probe')

    def prepare_subsequence(self):
        self.time_spinpol_mu = self.get_parameter('time_spinpol_us', group='timing', override=False, conversion=seconds_to_mu, units=us)
    
    @kernel(flags={"fast-math"})
    def run(self):
        # probe pulse
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
