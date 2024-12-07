from artiq.experiment import *
from LAX_exp.base import LAXDevice

import time
import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config


class Wavemeter(LAXDevice):
    """
    High-level API functions for configuring the wavemeter
    """
    name = "wavemeter"
    password = "lab"
    ip = "10.97.111.8"

    channels = {
        '397nm': (5,    '755.221845',   (0, 1), True,   5,  [-4, 4]),
        '423nm': (4,    '709.077640',   (0, 2), True,   4,  [-4, 4]),
        '729nm': (9,    '411.0416',     (0, 5), False,  -1, [-4, 4]),
        '854nm': (14,   '350.862460',   (0, 4), True,   8,  [-4, 4]),
        '866nm': (13,   '345.999945',   (0, 3), True,   7,  [-4, 4]),
    }

    alarm_threshold_mhz = 1e12

    def prepare_device(self):
        self.cxn = labrad.connect(self.ip, port=7682, tls_mode='off', username='', password='lab')
        self.wavemeter = self.cxn.multiplexerserver

    @rpc
    def read_channel_frequency(self, channel: TInt32) -> TFloat:
        """
        Read frequency of the wavemeter

        Args:
            channel (TInt32): wavemeter channel

        Returns
        """
        return self.wavemeter.get_frequency(channel)

    @rpc
    def set_channel_frequency(self, channel: TInt32, freq_thz: TFloat) -> TNone:
        """
        Set the frequency of the wavemeter

        Args:
            channel (TInt32): wavemeter channel
            freq_thz: frequency in THz to set wavemeter PID to
        """

        if channel == 5:
            assert 755.21 freq_thz < 755.23, "Check Set Frequency of 397 Channel for Wavemeter"


        self.wavemeter.set_pid_course(channel, freq_thz)