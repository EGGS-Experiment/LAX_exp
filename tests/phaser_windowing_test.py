from artiq.experiment import *
from artiq.coredevice.rtio import *
from artiq.coredevice.core import rtio_get_counter

import numpy as np
from numpy import int32, int64

_DMA_HANDLE_PHASER_TEST = "phaser_test_run"


class PhaserWindowingTest(EnvExperiment):
    """
    Phaser Windowing Test

    Test windowing on phaser0.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("phaser0")

        # general
        self.setattr_argument("repetitions",                            NumberValue(default=5, ndecimals=0, step=1, min=1, max=10000))
        
        # sequence timing
        self.setattr_argument("time_pulse_ms",                          NumberValue(default=5, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("time_sequence_delay_ms",                 NumberValue(default=5, ndecimals=5, step=1, min=1, max=10000))


    def prepare(self):
        # alias devices
        self.ph_channel =                                               self.phaser0.channel[0]

        # conversion values
        self.frequency_to_ftw =                                         (1 << 30) / 6.25e6
        self.pct_to_asf =                                               0x7FFF
        
        # timing values
        self.time_sample_mu =                                           int64(40)
        self.time_frame_mu =                                            self.phaser0.t_frame
        self.time_pulse_mu =                                            self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)
        self.time_sequence_delay_mu =                                   self.core.seconds_to_mu(self.time_sequence_delay_ms * ms)

        # ensure phaser pulse time is a multiple of the phaser frame period
        # (4 ns/clock * 8 clock cycles * 10 words = 320ns)
        if self.time_pulse_mu % self.time_frame_mu:
            t_frame_multiples = round(self.time_pulse_mu / self.time_frame_mu + 0.5)
            self.time_eggs_heating_mu = int64(self.time_frame_mu * t_frame_multiples)

        # set up windowing variables
        # double the sample period since we use two oscillators
        num_samples =                                                   self.time_eggs_heating_mu / (2 * self.time_sample_mu)
        max_amplitude =                                                 int32(0x7FFF / 2 - 1)

        # calculate windowing values - hann window
        self.ampl_window_mu_list =                                      np.power(np.sin(np.arange(num_samples) / num_samples * np.pi), 2)
        self.ampl_window_mu_list =                                      int32(self.ampl_window_mu_list * max_amplitude)


    @kernel(flags={'fast-math'})
    def run(self):
        # initialize
        self.core.reset()
        self.initialize_phaser()

        # record window
        self.window_record()
        handle_test = self.core_dma.get_handle(_DMA_HANDLE_PHASER_TEST)

        # continually output window
        for i in range(self.repetitions):
            
            # align to next phaser frame
            at_mu(self.phaser0.get_next_frame_mu())
            
            # output signal 
            self.core_dma.playback_handle(handle_test)
            
            # wait a holdoff period
            delay_mu(self.time_sequence_delay_mu)


    @kernel(flags={'fast-math'})
    def initialize_phaser(self):
        # initialize phaser board
        self.core.break_realtime()
        self.phaser0.init(debug=True)
        self.core.break_realtime()

        # trf setup, and disable rf output while we set things up
        self.ph_channel.set_att(0 * dB)
        self.ph_channel.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()

        # set nco to center eggs rf around 85 MHz exactly
        at_mu(self.phaser0.get_next_frame_mu())
        self.ph_channel.set_nco_frequency(-217.083495 * MHz)
        self.ph_channel.set_nco_phase(0.)
        self.phaser0.dac_sync()
        self.core.break_realtime()

        # duc
        self.ph_channel.set_duc_frequency(0 * MHz)
        self.ph_channel.set_duc_cfg()
        self.phaser0.duc_stb()
        self.core.break_realtime()

        # set attenuations
        self.ph_channel.set_att(0 * dB)
        self.core.break_realtime()

        # set oscillators (i.e. sidebands)
        self.ph_channel.oscillator[0].set_frequency(self.freq_sideband)
        self.ph_channel.oscillator[0].set_amplitude_phase(0., clr=0)
        self.core.break_realtime()
        self.ph_channel.oscillator[1].set_frequency(self.freq_eggs_heating_secular_mhz * MHz)
        self.ph_channel.oscillator[1].set_amplitude_phase(0., clr=0)
        self.core.break_realtime()

        # re-enable rf output
        self.ph_channel.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()

    @kernel(flags={'fast-math'})
    def window(self):
        # record windowing sequence
        with self.core_dma.record(_DMA_HANDLE_PHASER_TEST):
            
            # iterate across window amplitude values in mu
            for amp_val_mu in self.ampl_window_mu_list:

                # set amplitude for oscillator 0
                self.ph_channel.oscillator[0].set_amplitude_phase_mu(amp_val_mu)
                delay_mu(self.time_sample_mu)

                # set amplitude for oscillator 1
                self.ph_channel.oscillator[1].set_amplitude_phase_mu(amp_val_mu)
                delay_mu(self.time_sample_mu)

        self.core.break_realtime()

    def analyze(self):
        print('\t\tphaser test done')
