from artiq.experiment import *

import numpy as np
from datetime import datetime
from LAX_exp.extensions import *
from LAX_exp.base import LAXSequence


class ProbeAmplitudeCalibration(LAXSequence):
    """
    Sequence: Probe Amplitude Calibration

    Adjust the probe asf to meet a target value across a range of frequencies.
    """
    name = 'probe_amplitude_calibration'

    # todo: change everything from here on below
    parameters = {
        "sample_num":                       ('calibration.pmt.sample_num',              None),
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu),
        'time_readout_mu':                  ('timing.time_readout_us',                  us_to_mu),
    }
    devices = [
        'pump',
        'repump_cooling',
        'pmt'
    ]

    def prepare_sequence(self):
        # create subsequence datasets
        self.set_dataset("counts_signal", np.zeros(self.sample_num, dtype=np.int32))
        self.setattr_dataset("counts_signal")

        self.set_dataset("counts_noise", np.zeros(self.sample_num, dtype=np.int32))
        self.setattr_dataset("counts_noise")

    @kernel(flags={"fast-math"})
    def run(self):
        # set up 397 cooling
        self.pump.readout()
        self.core.break_realtime()

        # take discrimination counts w/cooling repump off
        self.pump.on()
        self.repump_cooling.off()
        self.core.break_realtime()

        for i in range(self.sample_num):
            self.pmt.count(self.time_readout_mu)
            self.counts_noise[i] = self.pmt.fetch_count()
            self.core.break_realtime()

        # take discrimination counts w/cooling repump on
        self.repump_cooling.on()
        self.core.break_realtime()

        for i in range(self.sample_num):
            self.pmt.count(self.time_readout_mu)
            self.counts_signal[i] = self.pmt.fetch_count()
            self.core.break_realtime()

    def analyze(self):
        # process count data
        signal_mean, noise_mean = (np.mean(self.counts_signal), np.mean(self.counts_noise))
        signal_std, noise_std = (np.std(self.counts_signal), np.std(self.counts_noise))

        # calculate discrimination threshold value
        threshold = signal_mean / np.log(1 + signal_mean / noise_mean)

        print('\tsignal mean: {} +/- {}'.format(signal_mean, signal_std))
        print('\tnoise mean: {} +/- {}'.format(noise_mean, noise_std))
        print('\tthreshold: {}'.format(threshold))

        # update values in scheduler
        self.set_dataset("calibration.pmt.counts_signal", signal_mean, broadcast=True, persist=True)
        self.set_dataset("calibration.pmt.counts_noise", noise_mean, broadcast=True, persist=True)
        self.set_dataset("calibration.pmt.calibration_timestamp", datetime.timestamp(datetime.now()), broadcast=True, persist=True)