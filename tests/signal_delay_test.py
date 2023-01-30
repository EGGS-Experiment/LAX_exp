from artiq.experiment import *


class SignalDelayTest(EnvExperiment):
    """
    Signal Delay Test

    Send out signals simultaneously to measure the delay.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # get devices
        self.setattr_argument("ttl_trigger_channel_num",                NumberValue(default=20, ndecimals=0, step=1, min=0, max=23))
        self.setattr_argument("dds_signal_board_num",                   NumberValue(default=0, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_signal_channel_num",                 NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

    def prepare(self):
        # get specified devices
        self.ttl_channel =                                              self.get_device('ttl{:d}'.format(self.ttl_trigger_channel_num))
        self.dds_channel =                                              self.get_device('urukul{:d}_ch{:d}'.format(self.dds_signal_board_num, self.dds_signal_channel_num))

    @kernel
    def run(self):
        self.core.reset()

        # prepare devices
        self.ttl_channel.off()
        self.dds_channel.cfg_sw(True)
        delay_mu(1000)

        # trigger signal
        self.ttl_channel.on()
        delay_mu(1000)

        # turn devices off again
        self.ttl20.off()
        self.urukul1_ch3.cfg_sw(False)

    def analyze(self):
        pass
