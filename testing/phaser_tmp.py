import labrad
import numpy as np
from os import environ
from artiq.experiment import *


class phasertesttmp(EnvExperiment):
    """
    phasertest tmp
    phaser test tmp
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device('phaser0')

    def prepare(self):
        pass

    @kernel(flags={'fast-math'})
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        delay_mu(1000000)
        
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(1000000)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.0, clr=1)
        delay_mu(1000000)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(0. * MHz)
        delay_mu(1000000)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg()
        delay_mu(1000000)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()
        delay_mu(1000000)

        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].oscillator[0].set_frequency(0. * MHz)
        delay_mu(1000000)



    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_nco_frequency(100. * MHz)
    #     # delay_mu(40)
    #     # self.phaser0.channel[1].set_nco_frequency(100. * MHz)
    #     # # clear DAC NCO phase
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_nco_phase(0.)
    #     # delay_mu(40)
    #     # self.phaser0.channel[1].set_nco_phase(0.)
    #     # # sync DAC for both channels
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.dac_sync()
    #     # # add slack
    #     # self.core.break_realtime()
    #
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].set_att(0. * dB)
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[1].set_att(0. * dB)
    #
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_duc_frequency(5. * MHz)
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[1].set_duc_frequency(5. * MHz)
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[0].set_duc_cfg()
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.channel[1].set_duc_cfg()
    #     # at_mu(self.phaser0.get_next_frame_mu())
    #     # self.phaser0.duc_stb()
    #
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].oscillator[0].set_frequency(5. * MHz)
    #     at_mu(self.phaser0.get_next_frame_mu())
    #     self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.999, phase=0., clr=0)
    #     self.core.break_realtime()
    #