from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad


class Wavemeter(LAXDevice):
    """
    High-level API functions for configuring the wavemeter
    """
    name = "wavemeter"
    kernel_invariants = {
        "cxn", "wavemeter",
        "password", "ip", "channels",
    }

    def prepare_device(self):
        # get wavemeter configs
        # todo: get from config file instead
        self.password = "lab"
        self.ip = "10.97.111.8"
        self.channels = {
            '397nm': (5, '755.221845', (0, 1), True, 5, [-4, 4]),
            '423nm': (4, '709.077640', (0, 2), True, 4, [-4, 4]),
            '729nm': (9, '411.0416', (0, 5), False, -1, [-4, 4]),
            '854nm': (14, '350.862460', (0, 4), True, 8, [-4, 4]),
            '866nm': (13, '345.999945', (0, 3), True, 7, [-4, 4]),
        }

        # create labrad connections
        self.cxn = labrad.connect(self.ip, port=7682, tls_mode='off', username='', password='lab')
        self.wavemeter = self.cxn.multiplexerserver

    @rpc
    def read_channel_frequency(self, channel: TInt32) -> TFloat:
        """
        Read frequency of the wavemeter
        :param channel: wavemeter channel
        :return: frequency of the wavemeter channel (in THz).
        """
        return self.wavemeter.get_frequency(channel)

    @rpc
    def set_channel_frequency(self, channel: TInt32, freq_thz: TFloat) -> TNone:
        """
        Set the frequency of the wavemeter
        :param channel: wavemeter channel to set.
        :param freq_thz: target wavemeter frequency setpoint (in THz).
        """
        pass
        # if channel == 5:
        #     assert 755.21 freq_thz < 755.23, "Check Set Frequency of 397 Channel for Wavemeter"

