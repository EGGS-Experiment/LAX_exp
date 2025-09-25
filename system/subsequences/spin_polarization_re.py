from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SpinPolarizationRE(LAXSubsequence):
    """
    Subsequence: Spin Polarization RE

    Optically pump the ion in the |S-1/2, m_j=-1/2> state and out of the |D=5/2, m_j=-5/2> using the \sigma_- 397 beam,
    866 beam, and 854 beam
    """
    name = 'Spin Polarization RE'
    kernel_invariants = {
        "time_spinpol_mu"
    }

    def build_subsequence(self):
        self.setattr_device('probe')
        self.setattr_device('repump_qubit')
        self.setattr_device('repump_cooling')

    def prepare_subsequence(self):
        self.time_spinpol_mu = self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                  conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # probe pulse
        self.probe.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        self.repump_qubit.off()
        self.repump_cooling.off()
