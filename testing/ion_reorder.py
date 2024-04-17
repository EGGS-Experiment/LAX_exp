import numpy as np
from artiq.experiment import *


class IonReorder(EnvExperiment):
    """
    Utility: Ion Reorder
    """
    kernel_invariants = {
        "dds", "dds_att", "dds_ampl", "dds_freq", "dds_profile", "dds_time", "dds_off_att"
    }


    def build(self):
        # get core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")

        self.setattr_device("ttl0_counter")
        self.setattr_device("ttl15")
        self.setattr_device("urukul1_ch3")


    def prepare(self):
        self.dds =          self.get_device("urukul1_ch3")
        self.dds_att =      self.dds.cpld.att_to_mu(0. * dB)
        self.dds_ampl =     self.dds.amplitude_to_asf(0.98)
        self.dds_freq =     self.dds.frequency_to_ftw(1375 * kHz)
        self.dds_profile =  4

        self.dds_time =     self.core.seconds_to_mu(1000 * ms)
        self.dds_off_att =  self.dds.cpld.att_to_mu(31.5 * dB)

    @kernel(flags={"fast-math"})
    def run(self):
        self.run_prepare()

        # shuffle
        self.shuffle()

        # clean up
        self.cleanup()

    @kernel(flags={"fast-math"})
    def run_prepare(self):
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.reset()

    @kernel(flags={"fast-math"})
    def shuffle(self):
        self.core.break_realtime()

        # set waveform
        self.dds.set_att_mu(self.dds_att)
        self.dds.set_mu(self.dds_freq, asf=self.dds_ampl, profile=self.dds_profile)
        self.core.break_realtime()

        # shuffle
        self.dds.cpld.set_profile(self.dds_profile)
        self.dds.sw.on()
        delay_mu(self.dds_time)
        self.dds.sw.off()

        # clean up to prevent DDS overheating
        self.dds.set_att_mu(self.dds_off_att)
        self.dds.set_mu(self.dds_freq, asf=0x01, profile=self.dds_profile)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup(self):
        self.core.break_realtime()

        # set urukul
        self.dds.sw.off()
        self.dds.set_att_mu(self.dds_off_att)
        self.dds.set_mu(self.dds_freq, asf=0x01, profile=self.dds_profile)

        # wait until end
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
