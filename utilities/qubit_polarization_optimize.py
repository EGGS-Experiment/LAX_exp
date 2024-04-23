import labrad
from os import environ
from artiq.experiment import *


class QubitPolarizationOptimize(EnvExperiment):
    """
    Utility: Qubit Polarization Optimize

    Optimize the polarization of the 729nm waveplate
    """


    def build(self):
        # todo: signal min, signal error
        # voltage values
        # self.dc_channel_dict =                      dc_config.channeldict
        # self.setattr_argument("dc_channel_name",    EnumerationValue(list(self.dc_channel_dict.keys()), default='V Shim'))
        self.setattr_argument("dc_voltage_v",       NumberValue(default=10.0, ndecimals=3, step=1, min=0, max=400))


    def prepare(self):
        # get voltage channel number
        # self.dc_channel_num =       self.dc_channel_dict[self.dc_channel_name]['num']

        # connect to labrad
        self.cxn =                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                   self.cxn.dc_server


    def run(self):
        # set DC voltage
        self.dc.voltage(self.dc_channel_num, self.dc_voltage_v)
