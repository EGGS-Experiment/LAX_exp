import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandReadout(LAXSubsequence):
    """
    Subsequence: Sideband Readout

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    Note: this subsequence handles profiles/values by itself.
    """
    name = 'sideband_readout'
    kernel_invariants = {
        "profile_dds",
        "ampl_sideband_readout_asf", "att_sideband_readout_mu", "time_sideband_readout_mu",
        "freq_sideband_readout_ftw_list"
    }

    def build_subsequence(self, profile_dds: TInt32 = 0) -> TNone:
        """
        Defines the main interface for the subsequence.
        Arguments:
            profile_dds: the AD9910 profile to use for sideband readout.
        """
        # set subsequence parameters
        self.profile_dds = profile_dds

        # sideband cooling readout
        self.setattr_argument("freq_rsb_readout_mhz_list",      Scannable(
                                                                    default=[
                                                                        CenterScan(100.7048, 0.03, 0.0006, randomize=True),
                                                                        ExplicitScan([100.7044]),
                                                                    ],
                                                                    global_min=30, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, precision=5
                                                                ), group=self.name)
        self.setattr_argument("freq_bsb_readout_mhz_list",      Scannable(
                                                                    default=[
                                                                        CenterScan(101.4081, 0.03, 0.0006, randomize=True),
                                                                        ExplicitScan([101.3901]),
                                                                    ],
                                                                    global_min=30, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, precision=5
                                                                ), group=self.name)
        self.setattr_argument("ampl_sideband_readout_pct",      NumberValue(default=50, precision=3, step=10, min=1, max=50., scale=1., unit="%"), group=self.name)
        self.setattr_argument("att_sideband_readout_db",        NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, scale=1., unit="dB"), group=self.name)
        self.setattr_argument("time_sideband_readout_us",       NumberValue(default=26.11, precision=5, step=1, min=1, max=10000, scale=1., unit="us"), group=self.name)

        # get relevant devices
        self.setattr_device('qubit')

    def prepare_subsequence(self):
        # prepare readout waveform values
        self.ampl_sideband_readout_asf =        self.qubit.amplitude_to_asf(self.ampl_sideband_readout_pct / 100.)
        self.att_sideband_readout_mu =          att_to_mu(self.att_sideband_readout_db * dB)
        self.time_sideband_readout_mu =         self.core.seconds_to_mu(self.time_sideband_readout_us * us)

        # combine readout frequencies WITHOUT shuffling them
        self.freq_sideband_readout_ftw_list =   np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                          for freq_mhz in (list(self.freq_rsb_readout_mhz_list) +
                                                                           list(self.freq_bsb_readout_mhz_list))])

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # set readout waveform for qubit
        self.qubit.set_profile(self.profile_dds)
        self.qubit.set_att_mu(self.att_sideband_readout_mu)

        # population transfer pulse
        self.qubit.on()
        delay_mu(self.time_sideband_readout_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def run_time(self, time_readout_mu: TInt64) -> TNone:
        # set readout waveform for qubit
        self.qubit.set_profile(self.profile_dds)
        self.qubit.set_att_mu(self.att_sideband_readout_mu)

        # population transfer pulse
        self.qubit.on()
        delay_mu(time_readout_mu)
        self.qubit.off()

