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

        self.setattr_argument("freq_nco_mhz", NumberValue(default=-217.083495, ndecimals=6, step=100, min=-400., max=400.))

    def prepare(self):
        # ensure NCO frequency is valid
        if (self.freq_nco_mhz > 400.) or (self.freq_nco_mhz < -400.):
            raise Exception("Invalid phaser NCO frequency. Must be in range [-400, 400].")
        elif (self.freq_nco_mhz > 300.) or (self.freq_nco_mhz < -300.):
            print("Warning: Phaser NCO frequency outside passband of [-300, 300] MHz.")

        # set relevant values for phaser initialization
        self.time_phaser_sample_mu = int64(40)


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        '''
        *************ATTENUATORS*******************
        '''
        # set maximum attenuation to eliminate output
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_att(31.5 * dB)

        # add slack
        self.core.break_realtime()


        '''
        *************PHASER*******************
        '''
        # initialize phaser
        self.phaser0.init(debug=True)

        # add slack
        self.core.break_realtime()


        '''
        *************DAC/NCO*******************
        '''
        # set DAC NCO frequency to center output at 85 MHz exactly
        # note: TRF372017 freq is 302.083918 MHz => DAC NCO should be 217.083918 MHz
        # note: currently using -217.083495 as NCO center? idk why
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_frequency(self.freq_nco_mhz * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_frequency(self.freq_nco_mhz * MHz)


        # clear DAC NCO phase
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_nco_phase(0.)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_nco_phase(0.)

        # sync DAC for both channels
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_sync()

        # add slack
        self.core.break_realtime()


        '''
        *************DUC (DIGITAL UPCONVERTER)*******************
        '''
        # set channel DUC frequencies to 0
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(0. * MHz)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_frequency(0. * MHz)

        # clear channel DUC phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)

        # strobe update register for both DUCs
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()

        # add slack
        self.core.break_realtime()


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

        # add slack
        self.core.break_realtime()


        '''
        *************TRF (UPCONVERTER)*******************
        '''
        # enable outputs for both channels here
        # note: want to do this at end instead of beginning since output may be nonzero and
        # will be cycling through frequencies as we initialize components
        # note: want to leave trf outputs persistently enabled since phase relation
        # between channels can change after adjusting the TRF
        # note: no need to set TRF frequency here since we already do this in device_db
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].en_trf_out(rf=1, lo=0)
        delay_mu(self.time_phaser_sample_mu)
        self.phaser0.channel[1].en_trf_out(rf=1, lo=0)

        # add slack
        self.core.break_realtime()

