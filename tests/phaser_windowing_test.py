from artiq.experiment import *
from artiq.coredevice.rtio import *
from artiq.coredevice.core import rtio_get_counter

import numpy as np
from numpy import int32, int64


class PhaserWindowingTest(EnvExperiment):
    """
    Phaser Windowing Test

    Test windowing on phaser0.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("phaser0")

    def prepare(self):
        # alias devices
        self.ph_channel = self.phaser0.channel[0]

        # base values
        self.t_sample_mu = int64(40)
        self.frequency_to_ftw = np.float64((1 << 30) / 6.25e6)
        self.pct_to_asf = 0x7FFF

        # configure windowing
        time_output_mu = int64(50000000)

        self.ampl_window_mu_list = np.array([], dtype=int32)

    @kernel(flags={'fast-math'})
    def run(self):
        self.core.reset()
        self.initialize_phaser()
        self.window()

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
        # iterate across window amplitudes
        for amp_val_mu in self.ampl_window_mu_list:
            self.ph_channel.oscillator[0].set_amplitude_phase_mu(amp_val_mu)
            delay_mu(self.t_sample_mu)
            self.ph_channel.oscillator[1].set_amplitude_phase_mu(amp_val_mu)
            delay_mu(self.t_sample_mu)
