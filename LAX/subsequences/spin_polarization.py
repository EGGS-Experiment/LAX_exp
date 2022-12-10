from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class SpinPolarization(LAXSubsequence):
    """
    Subsequence: Spin Polarization
        Place the ion in the S-1/2 mj=-1/2 state using the polarized 397 probe beam.
    """
    name = 'spin_polarization'

    parameters = {
        'time_spinpol_mu':                   ('timing.time_spinpol_us',                   us_to_mu)
    }
    devices = [
        'probe'
    ]
    
    @kernel(flags={"fast-math"})
    def run(self):
        # probe pulse
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
