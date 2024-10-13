import numpy as np
from artiq.experiment import *


class TicklePulse(EnvExperiment):
    """
    Utility: Tickle Pulse
    """
    kernel_invariants = {
        "dds",
        "dds_att_mu",
        "dds_freq_ftw", "dds_ampl_asf", "dds_profile",
        "dds_time_mu"
    }


    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("urukul1_ch3")

        # tickle arguments
        self.setattr_argument('freq_tickle_khz',    NumberValue(default=1372., precision=3, step=10, min=1, max=400000))
        self.setattr_argument('ampl_tickle_pct',    NumberValue(default=99., precision=2, step=10, min=0., max=100.))
        self.setattr_argument('att_tickle_db',      NumberValue(default=0., precision=1, step=0.5, min=0, max=31.5))
        self.setattr_argument('time_tickle_ms',     NumberValue(default=100., precision=3, step=10, min=0.01, max=100000))

    def prepare(self):
        # alias devices
        self.dds =          self.get_device("urukul1_ch3")

        # tickle waveform
        self.dds_att_mu =       self.dds.cpld.att_to_mu(self.att_tickle_db * dB)
        self.dds_freq_ftw =     self.dds.frequency_to_ftw(self.freq_tickle_khz * kHz)
        self.dds_ampl_asf =     self.dds.amplitude_to_asf(self.ampl_tickle_pct / 100.)
        self.dds_time_mu =      self.core.seconds_to_mu(self.time_tickle_ms * ms)
        self.dds_profile =      4


    @kernel(flags={"fast-math"})
    def run(self):
        # synchronize with timeline
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.reset()

        # shuffle & clean up
        self.shuffle()
        self.cleanup()

    @kernel(flags={"fast-math"})
    def shuffle(self):
        # set waveform
        self.dds.set_att_mu(self.dds_att_mu)
        self.dds.set_mu(self.dds_freq_ftw, asf=self.dds_ampl_asf, profile=self.dds_profile)
        self.core.break_realtime()

        # shuffle
        self.dds.cpld.set_profile(self.dds_profile)
        self.dds.sw.pulse_mu(self.dds_time_mu)

    @kernel(flags={"fast-math"})
    def cleanup(self):
        self.core.break_realtime()

        # clear urukul output to prevent leakage
        self.dds.sw.off()
        self.dds.set_att_mu(0)
        self.dds.set_mu(0x01, asf=0x01, profile=self.dds_profile)

        # wait until end
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

