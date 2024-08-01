import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

from LAX_exp.base import LAXExperiment, LAXEnvironment
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper
from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard

import matplotlib.pyplot as plt


class PulseShaperTest(LAXExperiment, Experiment):
    """
    PulseShaperTest
    """

    def build_experiment(self):
        # general
        self.setattr_argument("repetitions",            NumberValue(default=1, ndecimals=0, step=1, min=1, max=1000))

        # timing
        self.setattr_argument("time_reset_us",          NumberValue(default=2000, ndecimals=3, step=500, min=0.001, max=100000))

        # waveform
        self.setattr_argument("freq_carrier_mhz",       NumberValue(default=50., ndecimals=6, step=1., min=0., max=200.))
        self.setattr_argument("freq_sideband_khz",      NumberValue(default=1400, ndecimals=3, step=100, min=-10000, max=10000))
        self.setattr_argument("phase_ch1_turns",        NumberValue(default=0., ndecimals=3, step=0.1, min=-1.0, max=1.0))
        self.setattr_argument("att_phaser_db",          NumberValue(default=5., ndecimals=1, step=0.5, min=0, max=31.5))

        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")
        self.setattr_device('phaser_eggs')

        # objects
        self.spinecho_wizard = SpinEchoWizard(self)
        self.pulse_shaper = PhaserPulseShaper(self)

    def prepare_experiment(self):
        # prepare hardware values
        self.time_reset_mu =    self.core.seconds_to_mu(self.time_reset_us * us)
        self.freq_carrier_hz =  self.freq_carrier_mhz * MHz
        self.freq_sideband_hz = self.freq_sideband_khz * kHz
        # prepare waveform playback
        self._wav_idx = 0

        # create waveform
        self.spinecho_wizard.prepare()
        self.spinecho_wizard.calculate_pulseshape()
        self.spinecho_wizard.compile_waveform()
        # get waveform data
        self._wav_data_ampl, self._wav_data_phas, self._wav_data_time = self.spinecho_wizard.get_waveform()

    @property
    def results_shape(self):
        return (2, 2)

    @kernel(flags={"fast-math"})
    def run_main(self):
        # setup hardware
        self.run_prepare()

        # run waveform
        for i in range(self.repetitions):
            # play waveform
            at_mu(self.phaser_eggs.get_next_frame_mu())
            self.ttl8.on()
            self.pulse_shaper.waveform_playback(self._wav_idx)
            self.ttl8.off()

            # delay reset time
            delay_mu(self.time_reset_mu)

    @kernel(flags={"fast-math"})
    def run_prepare(self):
        self.core.break_realtime()
        self.core.reset()

        # prepare TTLs
        self.ttl8.off()
        self.ttl9.off()

        # prepare phaser - attenuators
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(self.att_phaser_db * dB)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att(self.att_phaser_db * dB)
        self.core.break_realtime()

        # prepare phaser - frequencies
        self.phaser_configure(self.freq_carrier_hz, self.freq_sideband_hz)
        self.core.break_realtime()

        # waveform - record
        # delay_mu(1000000)
        self._wav_idx = self.pulse_shaper.waveform_record(self._wav_data_ampl, self._wav_data_phas, self._wav_data_time)
        self.core.break_realtime()

        # waveform - load
        self.pulse_shaper.waveform_load()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat) -> TNone:
        """
        todo: document
        """
        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(self.phase_ch1_turns)
        self.phaser_eggs.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)

    def analyze(self):
        self.spinecho_wizard.display_waveform()
