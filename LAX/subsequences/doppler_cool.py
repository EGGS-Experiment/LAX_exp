import numpy as np
from artiq.experiment import *
# todo: allow values to be specified, but if not, take from dataset
# todo: get relevant devices
# todo: create function that
# todo: need sequence base class as well
# core DMA should only be for sequence level, otherwise overhead of recording will likely pile up


class DopplerCooling(HasEnvironment):
    """
    Sequence: Doppler Cooling
        Apply Doppler cooling for a given time.
    """

    kernel_invariants = {}


    def build(self, cooling_time_us=None):
        """
        Set devices and arguments for the experiment.
        """
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # get parameters
        if cooling_time_us is None:
            cooling_time_us = self.get_dataset()

        setattr(
            self,
            "time_doppler_cooling_mu",
            self.core.seconds_to_mu(cooling_time_us * us)
        )


    @kernel(flags={"fast-math"})
    def run(self, cooling_time_mu=self.time_doppler_cooling_mu):
        """
        MAIN SEQUENCE.

        Set the main Urukul board to cooling (profile 0), then pulse 397nm pump for given time.

        Arguments:
            cooling_time_mu     (np.int64)  : the time to doppler cool for, in machine units.
        """
        # set cooling waveform
        with parallel:
            self.dds_board.set_profile(0)
            delay_mu(self.time_profileswitch_delay_mu)

        # doppler cooling
        self.dds_board.cfg_switches(0b0110)
        delay_mu(self.time_doppler_cooling_mu)
        self.dds_board.cfg_switches(0b0100)
