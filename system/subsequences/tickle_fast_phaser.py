import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class TickleFastPhaser(LAXSubsequence):
    """
    Subsequence: Tickle Fast Phaser

    Heat the ion by applying an RF signal from a DDS at the secular frequency.
    Interferes two channels destructively to achieve fast switching times.
    """
    name = 'tickle_fast_phaser'


    def build_subsequence(self):
        self.setattr_argument('att_ticklefast_phaser_db', NumberValue(default=5, ndecimals=1, step=0.5, min=0, max=31.5), group='ticklefast')

        # get relevant devices
        self.setattr_device('phaser_eggs')

        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        # tmp remove

    def prepare_subsequence(self):
        # get DDS configuration parameters from dataset manager
        self.ampl_ticklefast_pct =          self.get_parameter('ampl_ticklefast_pct',
                                                               group='dds.ampl_pct', override=False)

        # prepare parameters for tickle pulse
        self.ampl_ticklefast_asf =          self.phaser_eggs.amplitude_to_asf(self.ampl_ticklefast_pct / 100.)
        self.att_ticklefast_mu =            att_to_mu(self.att_ticklefast_phaser_db * dB)

        # set empty holder variable for delay time
        self.time_delay_mu =                np.int64(0)

        # set empty variables for phase
        self.phase_final_pow =              np.int32(0)

        # tmp remove
        self.time_system_prepare_delay_mu = np.int64(1850)
        # tmp remove

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # get starting time
        time_start_mu = self.phaser_eggs.get_next_frame_mu()

        # set phaser attenuation
        at_mu(time_start_mu)
        self.phaser_eggs.channel[0].set_att_mu(self.att_ticklefast_mu)

        # ensure DUC phase offset is cleared
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].set_duc_phase_mu(0x0)

        # set oscillator 0 frequency to 0 Hz (since frequency hopping will be done by DUC)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[0].oscillator[0].set_frequency_mu(0x0)

        # strobe DUC update register
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.duc_stb()

        # tmp remove
        self.ttl8.off()
        # tmp remove


    @kernel(flags={"fast-math"})
    def run(self):
        # get fiduciary start time
        # time_start_mu = self.phaser_eggs.get_next_frame_mu()

        # todo: do we need to clear duc phase?

        # start phaser output
        at_mu(self.phaser_eggs.get_next_frame_mu())
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase_mu(asf=self.ampl_ticklefast_asf, pow=self.phase_final_pow, clr=0)

            # tmp remove
            with sequential:
                delay_mu(self.time_system_prepare_delay_mu)
                self.ttl8.on()
                delay_mu(40)
                self.ttl8.off()
            # tmp remove

            # set time to next phaser sample period
            delay_mu(self.phaser_eggs.t_sample_mu)

        # stop phaser output
        self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase_mu(asf=0x0, pow=0x0, clr=1)


    @kernel
    def configure(self, freq_ftw: TInt32, phase_pow: TInt32, time_delay_mu: TInt64):
        # ensure timing is a multiple of 40ns (phaser sample period)
        self.time_delay_mu =                time_delay_mu
        if time_delay_mu % self.phaser_eggs.t_sample_mu:
            # round eggs heating time up to the nearest multiple of phaser frame period
            t_sample_multiples =            round(self.time_delay_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_delay_mu =            np.int64(self.phaser_eggs.t_sample_mu * t_sample_multiples)

        # store phase offset word
        self.phase_final_pow =              phase_pow
        self.core.break_realtime()

        # set waveforms for phaser and strobe update register
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency_mu(freq_ftw - self.phaser_eggs.freq_center_ftw)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.duc_stb()

    def analyze(self):
        pass
        # print('\n\tresults:')
        # print('\t\turukul0 profile 1 phase: {:.4f}'.format(self.dds_ch0.pow_to_turns(self.phase_ch1_final_pow)))
        # print('\n')
