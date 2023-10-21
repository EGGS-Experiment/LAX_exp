from numpy import int64
from artiq.experiment import *


class PhaserInitialize(EnvExperiment):
    """
    Utility: Phaser Initialize

    Initialize the phaser.
    """

    def build(self):
        # set devices for initialization
        self.setattr_device("core")
        self.setattr_device("phaser0")

    def prepare(self):
        # set relevant values for phaser initialization
        self.time_phaser_sample_mu = int64(40)

    @kernel
    def run(self):
        self.core.reset()

        '''
        *************PHASER*******************
        '''
        # initialize phaser
        self.phaser0.init(debug=True)


        '''
        *************DAC*******************
        '''
        # set DAC NCO frequency to center output at 85 MHz exactly
        # todo note: trf frequency is 302.083918
        # todo note: 217.083918 is the right frequency
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_frequency((-217.083495) * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_frequency((-217.083495) * MHz)


        # clear DAC NCO phase
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_phase(0.)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_phase(0.)

        # sync DAC for both channels
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()


        '''
        *************DUC (DIGITAL UPCONVERTERS)*******************
        '''
        # set channel DUC frequencies to 0
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(0 * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_frequency(0 * MHz)

        # clear channel DUC phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)

        # strobe update register for both DUCs
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()


        '''
        *************OSCILLATORS*******************
        '''
        # reset oscillator frequency and amplitude, and keep phase accumulator persistently cleared
        # note: this has to happen before TRF or attenuator adjustment to ensure channel outputs are 0
        for i in range(5):
            # clear channel 0 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[0].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)

            # clear channel 1 oscillator
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[1].oscillator[i].set_frequency(0 * MHz)
            delay_mu(self.time_phaser_sample_mu)
            self.phaser0.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., clr=1)


        '''
        *************TRF (UPCONVERTER)*******************
        '''
        # enable outputs for both channels
        # note: want to leave trf outputs persistently enabled since phase relation
        # between channels can change after adjusting the TRF
        # note: no need to set TRF frequency here since we already do this in device_db
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].en_trf_out(rf=1, lo=0)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].en_trf_out(rf=1, lo=0)


        '''
        *************ATTENUATORS*******************
        '''
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_att(31.5 * dB)

        self.core.break_realtime()
