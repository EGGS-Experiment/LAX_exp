from artiq.experiment import *

import numpy as np
from datetime import datetime
from LAX_exp.extensions import *
from LAX_exp.base import LAXSequence


class PMTDiscrimination(LAXSequence):
    """
    Sequence: PMT Discrimination

    Record incoming photon counts to discriminate between the 0 and 1 states.
    """
    name = 'pmt_discrimination'

    parameters = {
        "calibration.pmt.sample_num":       ('calibration.pmt.sample_num',              100),
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu),
        'time_readout_mu':                  ('timing.time_readout_us',                  us_to_mu),
    }
    devices = [
        'pump',
        'repump_cooling',
        'pmt'
    ]

    # sequence parameters
    counts_signal =             []
    counts_noise =              []

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.break_realtime()

        # set up cooling
        self.pump.cooling()

        # take discrimination counts w/cooling repump off
        self.repump_cooling.off()

        # read out on PMT
        for i in range(self.num_reads):
            self.counts_noise.append(self.pmt.count(self.time_readout_mu))
            self.core.break_realtime()

        # take discrimination counts w/cooling repump on
        self.repump_cooling.on()

        # read out on PMT
        for i in range(self.num_reads):
            self.counts_signal.append(self.pmt.count(self.time_readout_mu))
            self.core.break_realtime()

    def analyze(self):
        # process count data
        signal_mean, noise_mean = (np.mean(self.counts_signal), np.mean(self.counts_noise))
        signal_std, noise_std = (np.std(self.counts_signal), np.std(self.counts_noise))

        # calculate discrimination threshold value
        threshold = signal_mean / np.log(1 + signal_mean / noise_mean)

        # update values in scheduler
        self.set_dataset("calibration.pmt.counts_signal", signal_mean, broadcast=True, persist=True)
        self.set_dataset("calibration.pmt.counts_noise", noise_mean, broadcast=True, persist=True)
        self.set_dataset("calibration.pmt.calibration_timestamp", datetime.timestamp(datetime.now()), broadcast=True, persist=True)
