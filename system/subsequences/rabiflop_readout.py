import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class RabiflopReadout(LAXSubsequence):
    """
    Subsequence: Rabi Flop Readout

    Read out by rabi flopping at a given transition.
    """
    name = 'rabiflop_readout'
    kernel_invariants = {
        "time_readout_mu_list",
        "freq_rabiflop_readout_ftw", "att_rabiflop_readout_mu", "ampl_qubit_asf"
    }

    def build_subsequence(self):
        # timing
        self.setattr_argument("time_readout_us_list",       Scannable(
                                                                default=RangeScan(0, 50, 51, randomize=True),
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group=self.name)
        # readout waveform parameters
        self.setattr_argument("freq_rabiflop_readout_mhz",  NumberValue(default=103.3455, precision=5, step=1, min=1, max=10000), group=self.name)
        self.setattr_argument("att_readout_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5), group=self.name)
        self.setattr_device('qubit')

    def prepare_subsequence(self):
        self.time_readout_mu_list =         np.array([self.core.mu_to_seconds(time_us * us)
                                                      for time_us in self.time_readout_us_list])

        self.freq_rabiflop_readout_ftw =    self.qubit.frequency_to_ftw(self.freq_rabiflop_readout_mhz * MHz)
        self.att_rabiflop_readout_mu =      att_to_mu(self.att_readout_db * dB)

        # get DDS amplitude for rabi flopping from dataset manager
        self.ampl_qubit_asf =   self.get_parameter('ampl_qubit_pct',
                                                   group='beams.ampl_pct', override=True,
                                                   conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        # set up the rabiflop readout profile
        self.qubit.set_mu(self.freq_rabiflop_readout_ftw, asf=self.ampl_qubit_asf,
                          profile=0, phase_mode=PHASE_MODE_CONTINUOUS)

    @kernel(flags={"fast-math"})
    def run(self, time_rabiflop_mu: TInt64) -> TNone:
        # set readout waveform for qubit
        self.qubit.set_profile(0)
        self.qubit.set_att_mu(self.att_rabiflop_readout_mu)

        # population transfer pulse
        self.qubit.on()
        delay_mu(time_rabiflop_mu)
        self.qubit.off()
