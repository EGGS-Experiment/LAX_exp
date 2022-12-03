from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu


class DopplerCool(LAXSubsequence):
    """
    Subsequence: Doppler Cooling
        Apply the 397nm pump beam for a given time.
    """
    name = 'doppler_cool'

    devices = [
        'pump'
    ]
    subsequence_parameters = {
        'time_doppler_cooling_mu':      ('timing.time_doppler_cooling_us', us_to_mu),
        'time_profileswitch_delay_mu':  ('timing.time_profileswitch_delay_us', us_to_mu)
    }

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveform
        with parallel:
            self.dds_board.set_profile(0)
            delay_mu(self.time_profileswitch_delay_mu)

        # doppler cooling
        self.pump.cfg_sw(1)
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.cfg_sw(0)
