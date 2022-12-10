from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class CorrelatePhotons(LAXSequence):
    """
    Sequence: Correlate Photons
        Correlate incoming photons with a synchronization signal.
    """
    name = 'correlate_photons'

    parameters = {
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu)
    }
    devices = [
        'pump'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        self.pump.cooling()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()
