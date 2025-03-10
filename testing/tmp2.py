import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testarg342(EnvExperiment):
    """
    testarg342
    Testing
    """
    kernel_invariants = {
        "phaser"
    }

    def build(self):
        self.setattr_device("core")
        self.setattr_device("phaser1")

    def prepare(self):
        self.phaser = self.get_device('phaser1')

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset
        self.core.reset()

        # configure atts
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_att(31.5)

        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(0. * MHz)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_nco_phase(0.)

        # configure freq
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(101.05 * MHz)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_cfg()
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.duc_stb()


        # configure ampl & freq
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].oscillator[0].set_frequency(0. * MHz)

        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].oscillator[0].set_amplitude_phase(0.0)
