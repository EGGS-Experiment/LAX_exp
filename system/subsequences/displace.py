import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class Displace(LAXSubsequence):
    """
    Subsequence: Displace

    Displace the ion into a coherent state by applying a dipole tone
    close to the secular frequency.
    Differs from the other "tickle_<x>" subsequences since it is designed to be integrated with other
    motional subsequences (e.g. squeezing) and so ensures the hardware is always correctly set.
    """
    name = 'displace'
    kernel_invariants = {
        "att_displace_mu",
        "displace_func"
    }

    def build_subsequence(self):
        self.setattr_argument("enable_displacement",        BooleanValue(default=True), group=self.name)
        self.setattr_argument('att_displace_db',            NumberValue(default=10., precision=1, step=0.5, min=0, max=31.5), group=self.name)

        # get relevant devices
        self.setattr_device('dds_dipole')

    def prepare_subsequence(self):
        # prepare parameters for tickle pulse
        self.att_displace_mu =                              att_to_mu(self.att_displace_db * dB)

        # create empty holder variables to support later configuration
        self.freq_displace_ftw =                            np.int32(0)
        self.phase_displace_pow =                           np.int32(0)
        self.time_displace_mu =                             np.int64(0)

        # configure displacement on/off
        self.displace_func =                                self.dds_dipole.on
        if not self.enable_displacement:
            self.displace_func =                            self.dds_dipole.off

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set dds attenuation here - ensures that dds channel will have correct attenuation
        # for any sequences recorded into DMA during initialize_experiment
        self.dds_dipole.set_att_mu(self.att_displace_mu)

    @kernel(flags={"fast-math"})
    def configure(self, freq_ftw: TInt32, phase_pow: TInt32, time_mu: TInt64):
        # store parameters as instance attribute
        self.time_displace_mu = time_mu
        print(self.time_displace_mu)
        self.core.break_realtime()

        # set waveforms for profiles
        # displacement waveform
        self.dds_dipole.set_mu(freq_ftw, asf=self.dds_dipole.ampl_modulation_asf, profile=0,
                               pow_=phase_pow, phase_mode=PHASE_MODE_CONTINUOUS)
        # blank waveform
        self.dds_dipole.set_mu(freq_ftw, asf=0x0, profile=2,
                               pow_=phase_pow, phase_mode=PHASE_MODE_CONTINUOUS)

    @kernel(flags={"fast-math"})
    def run(self):
        # set blank output waveform
        self.dds_dipole.set_profile(2)
        self.dds_dipole.write32(_AD9910_REG_CFR1,
                                (1 << 16) |  # select_sine_output
                                (1 << 13))  # phase_autoclear

        # align to coarse RTIO clock
        time_start_mu = now_mu() & ~0x7
        # begin output waveform
        at_mu(time_start_mu)
        self.dds_dipole.set_profile(0)

        # open RF switches early since they have a ~100 ns rise time
        at_mu(time_start_mu
              + (TIME_URUKUL_BUS_WRITE_DELAY_MU + TIME_AD9910_PROFILE_SWITCH_DELAY_MU)
              - TIME_URUKUL_RFSWITCH_DELAY_MU)
        self.displace_func()

        # displace for given time
        # # note: turn off switch early to account for switch rise time
        # delay_mu(self.time_displace_mu
        #          - TIME_URUKUL_RFSWITCH_DELAY_MU)
        delay_mu(self.time_displace_mu)
        self.dds_dipole.off()
