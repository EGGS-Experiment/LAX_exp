from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

# tmp remove
import numpy as np
# tmp remove


class AbsorptionProbe(LAXSubsequence):
    """
    Subsequence: Absorption Probe

    Cool the ion before switching to a short, low-power pulse to measure absorption on a transition.
    """
    name = 'absorption_probe'

    def build_subsequence(self):
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

        # tmp remove
        self.setattr_device('sampler0')

    def prepare_subsequence(self):
        self.time_doppler_cooling_mu =      self.get_parameter('time_doppler_cooling_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_probe_mu =                self.get_parameter('time_probe_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # tmp remove
        self.reading_avg_tmp = 0
        self.set_dataset('sampler_tmp', [])
        # tmp remove

    # tmp remove
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        self.sampler0.set_gain_mu(6, 3)
    # tmp remove

    @kernel(flags={"fast-math"})
    def run(self):
        # tmp remove
        buffer_sampler = [0] * 8
        # tmp remove

        # set cooling waveform and activate cooling repump
        self.pump.cooling()
        self.repump_cooling.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # switch to probe waveform
        self.pump.set_profile(1)

        # record fluorescence of probe beam
        # self.pump.on()
        # self.pmt.count(self.time_probe_mu)
        # self.pump.off()

        # tmp remove
        with parallel:
            with sequential:
                self.pump.on()
                self.pmt.count(self.time_probe_mu)
                self.pump.off()
            with sequential:
                delay_mu(np.int64(self.time_probe_mu/2))
                self.sampler0.sample_mu(buffer_sampler)
                self.reading_avg_tmp = buffer_sampler[6]
                # read0 = buffer_sampler[6]
        # tmp remove

        # set cooling waveform and disable cooling repump
        self.pump.cooling()
        self.repump_cooling.off()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # switch to probe waveform
        self.pump.set_profile(1)

        # record fluorescence of probe beam
        # self.pump.on()
        # self.pmt.count(self.time_probe_mu)
        # self.pump.off()

        # tmp remove
        with parallel:
            with sequential:
                self.pump.on()
                self.pmt.count(self.time_probe_mu)
                self.pump.off()
            with sequential:
                delay_mu(np.int64(self.time_probe_mu/2))
                self.sampler0.sample_mu(buffer_sampler)
                self.reading_avg_tmp += buffer_sampler[6]
                self.append_to_dataset('sampler_tmp', self.reading_avg_tmp >> 1)
