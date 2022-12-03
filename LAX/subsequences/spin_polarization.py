from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class SpinPolarization(LAXSubsequence):
    """
    Subsequence: Spin Polarization
        Apply the 397nm probe beam for a given time.
    """
    name = 'spin_polarization'

    devices = [
        'probe'
    ]
    subsequence_parameters = {
        'time_redist_mu':           ('timing.time_redist_us', us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # readout pulse
        self.probe.cfg_sw(1)
        delay_mu(self.time_redist_mu)
        self.probe.cfg_sw(0)
