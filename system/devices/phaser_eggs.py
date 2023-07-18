from numpy import int32, int64
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
# todo: prepare - sample time, frame time
# todo: on/off func
# todo: att func


class PhaserEGGS(LAXDevice):
    """
    Device: Phaser (EGGS)

    Use the Phaser AWG to create the EGGS RF
    """
    name =          "phaser_eggs"
    core_device =   ('phaser', 'phaser0')

    def build_device(self):
        # alias both phaser output channels
        self.channel =                      self.phaser0.channel

        # set phaser sample/frame timings
        self.t_sample_mu =                  int64(40)
        self.t_frame_mu =                   int64(320)

    def prepare_device(self):
        # get frequency parameters
        # todo: get carrier center freq from dashboard instead
        self.freq_center_mhz =              85. * MHz
        self.freq_center_ftw =              int32(hz_to_ftw(self.freq_center_mhz * MHz))

        # get phase delay parameters
        # todo: inherent ch1 phase delay
        # todo: system ch1 phase delay
        self.phase_inherent_ch1_turns =     self.get_parameter('phas_ch1_inherent_turns', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)
        self.phase_system_ch1_turns =       self.get_parameter('phas_ch1_system_turns', group='beams.freq_mhz', override=False, conversion_function=hz_to_ftw, units=MHz)

    @kernel(flags={'fast-math'})
    def initialize_device(self):
        """
        Clear the DUC phase accumulators and sync the DAC.
        """
        # clear channel DUC phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.t_sample_mu)
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)

        # strobe update register for both DUCs
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()

        # sync DAC for both channels
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()

        # todo: set carrier frequency via DAC NCO frequency for both channels


    @kernel(flags={"fast-math"})
    def reset_oscillators(self):
        """
        Reset frequency and amplitude for all oscillators on both channels of the phaser.
        """
        for i in range(5):
            # clear channel 0 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[0].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[0].oscillator[i].set_amplitude_phase(amplitude=0.)

            # clear channel 1 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[1].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[1].oscillator[i].set_amplitude_phase(amplitude=0.)
