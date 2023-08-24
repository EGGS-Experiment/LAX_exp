import numpy as np
from artiq.experiment import *


class FastinoTMP1(EnvExperiment):
    """
    Fastino TMP1
    Set a value on the Fastino.
    """

    def build(self):
        # get devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("fastino0")

        # arguments
        self.setattr_argument("channel", NumberValue(default=0, ndecimals=0, step=1, min=0, max=32))
        self.setattr_argument("voltage", NumberValue(default=0, ndecimals=3, step=1, min=-10, max=10))

    def prepare(self):
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        self.setattr_device('phaser0')

        self.time_pulse_ms =            1000
        self.num_samples =              200
        self.time_rolloff_us =          100

        self.time_pulse_mu =            self.core.seconds_to_mu(self.time_pulse_ms * ms)
        self.time_rolloff_mu =          self.core.seconds_to_mu(self.time_rolloff_us * us)
        self.t_sample_mu =              np.int64(self.time_rolloff_mu / self.num_samples)

        self.ampl_arr =                 np.power(np.sin((np.pi / (2. * self.num_samples)) * np.linspace(1, self.num_samples, self.num_samples)), 2)
        self.ampl_rev_arr =             self.ampl_arr[:: -1]

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # set up devices
        with parallel:
            self.ttl8.off()
            self.ttl9.off()
            self.fastino0.set_continuous(0x00)
        self.core.break_realtime()

        # record and retrieve DMA sequence
        self.record_dma()
        handle_rise_dma = self.core_dma.get_handle('FASTINO_TMP1')
        handle_fall_dma = self.core_dma.get_handle('FASTINO_TMP2')
        self.core.break_realtime()

        # main
        self.core.break_realtime()
        with parallel:
            with sequential:
                delay_mu(2700)
                self.ttl8.on()
            self.core_dma.playback_handle(handle_rise_dma)

        delay_mu(self.time_pulse_mu)
        self.core_dma.playback_handle(handle_fall_dma)


        # cleanup
        self.core.break_realtime()
        self.fastino0.set_dac(self.channel, 0.)
        with parallel:
            self.ttl8.off()
            self.ttl9.off()

    @kernel(flags={"fast-math"})
    def record_dma(self):
        self.core.break_realtime()

        with self.core_dma.record('FASTINO_TMP1'):
            for volt_v in self.ampl_arr:
                self.fastino0.set_dac(self.channel, volt_v)
                delay_mu(self.t_sample_mu)

        self.core.break_realtime()
        with self.core_dma.record('FASTINO_TMP2'):
            for volt_v in self.ampl_rev_arr:
                self.fastino0.set_dac(self.channel, volt_v)
                delay_mu(self.t_sample_mu)

    def analyze(self):
        pass
