import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandReadout(LAXSubsequence):
    """
    Subsequence: Sideband Readout

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    """
    name = 'sideband_readout'
    kernel_invariants = {
        "ampl_sideband_readout_asf",
        "att_sideband_readout_mu",
        "time_sideband_readout_mu",
        "freq_sideband_readout_ftw_list"
    }

    def build_subsequence(self):
        # sideband cooling readout
        self.setattr_argument("freq_rsb_readout_mhz_list",              Scannable(
                                                                            default=[
                                                                                ExplicitScan([102.5025]),
                                                                                CenterScan(102.5025, 0.010, 0.0005, randomize=True)
                                                                            ],
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group=self.name)
        self.setattr_argument("freq_bsb_readout_mhz_list",              Scannable(
                                                                            default=[
                                                                                ExplicitScan([103.2610]),
                                                                                CenterScan(103.2610, 0.010, 0.0005, randomize=True)
                                                                            ],
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group=self.name)
        self.setattr_argument("ampl_sideband_readout_pct",              NumberValue(default=50, ndecimals=3, step=10, min=1, max=50.), group=self.name)
        self.setattr_argument("att_sideband_readout_db",                NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)
        self.setattr_argument("time_sideband_readout_us",               NumberValue(default=70, ndecimals=5, step=1, min=1, max=10000), group=self.name)
        self.setattr_device('qubit')

    def prepare_subsequence(self):
        self.ampl_sideband_readout_asf =                                self.qubit.amplitude_to_asf(self.ampl_sideband_readout_pct / 100.)
        self.att_sideband_readout_mu =                                  att_to_mu(self.att_sideband_readout_db * dB)
        self.time_sideband_readout_mu =                                 self.core.seconds_to_mu(self.time_sideband_readout_us * us)

        # combine & shuffle readout frequencies
        self.freq_sideband_readout_ftw_list =                           np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                                        for freq_mhz in (list(self.freq_rsb_readout_mhz_list) + list(self.freq_bsb_readout_mhz_list))])
        np.random.shuffle(self.freq_sideband_readout_ftw_list)

    @kernel(flags={"fast-math"})
    def run(self):
        # set readout waveform for qubit
        self.qubit.set_profile(0)
        self.qubit.set_att_mu(self.att_sideband_readout_mu)

        # population transfer pulse
        self.qubit.on()
        delay_mu(self.time_sideband_readout_mu)
        self.qubit.off()
