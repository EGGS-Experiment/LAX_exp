from artiq.experiment import *
from LAX_exp.base import LAXDevice
from EGGS_labrad.config.multiplexerclient_config import multiplexer_config

import labrad


class Wavemeter(LAXDevice):
    """
    High-level API functions for configuring the wavemeter
    """
    name = "wavemeter"
    password = "lab"
    ip = "10.97.111.8"

    channels = multiplexer_config.channels
    alarm_threshold_mhz = 1e12

    def prepare_device(self):
        self.cxn = labrad.connect(self.ip, port=7682, tls_mode='off', username='', password='lab')
        self.wavemeter = self.cxn.multiplexerserver

    @rpc
    def get_channel_frequency(self, channel: TInt32) -> TFloat:
        """
        Read frequency of the wavemeter

        Args:
            channel (TInt32): wavemeter channel

        Returns
        """
        return self.wavemeter.get_frequency(channel)

    @rpc
    def get_channel_lock_frequency(self, channel: TInt32) -> TFloat:
        """
        Read lock frequency of the wavemeter

        Args:
            channel (TInt32): wavemeter channel

        Returns
        """
        return self.wavemeter.get_pid_course(channel)

    @rpc
    def set_channel_lock_frequency(self, channel: TInt32, freq_thz: TFloat) -> TNone:
        """
        Set the frequency of the wavemeter

        Args:
            channel (TInt32): wavemeter channel
            freq_thz: frequency in THz to set wavemeter PID to
        """
        pass

        # if channel == 5:
        #     assert 755.21 freq_thz < 755.23, "Check Set Frequency of 397 Channel for Wavemeter"

