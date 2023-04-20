import labrad
import numpy as np
from time import sleep

from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from EGGS_labrad.config.dc_config import dc_config


class DCSet(EnvExperiment):
    """
    DC Set

    Set the
    """
    kernel_invariants = {
        'time_pmt_gating_mu',
        'dc_micromotion_channel',
        'dc_micromotion_voltage_v'
    }

    global_parameters = []


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_channel",                         EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_voltage_v",                       NumberValue(default=37.0, ndecimals=3, step=1, min=0, max=400))


    def prepare(self):
        # get voltage parameters
        self.dc_micromotion_channel_num =                           self.dc_micromotion_channeldict[self.dc_micromotion_channel]['num']
        self.dc_micromotion_channel_name =                          self.dc_micromotion_channel
        self.dc_micromotion_voltages_v_list =                       np.array(list(self.dc_micromotion_voltages_v_list))

        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                   self.cxn.dc_server


    def run(self):
        # set DC voltage
        self.voltage_set(self.dc_micromotion_channel_num, voltage_v)
        self.core.break_realtime()
