import labrad
from os import environ
from artiq.experiment import *
from EGGS_labrad.config.dc_config import dc_config


class DCSet(EnvExperiment):
    """
    Tool: DC Set

    Set a voltage on the AMO8 HV DC box.
    For use when chaining a series of experiments.
    """

    def build(self):
        # voltage values
        self.dc_channel_dict = dc_config.channeldict
        self.setattr_argument("dc_channel_name",    EnumerationValue(list(self.dc_channel_dict.keys()), default='V Shim'))
        self.setattr_argument("dc_voltage_v",       NumberValue(default=10.0, precision=3, step=1, min=0, max=400, scale=1., unit="V"))

    def prepare(self):
        # get voltage channel number
        self.dc_channel_num = self.dc_channel_dict[self.dc_channel_name]['num']

        # connect to labrad
        self.cxn =  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =   self.cxn.dc_server

    def run(self):
        # set DC voltage
        self.dc.voltage(self.dc_channel_num, self.dc_voltage_v)

