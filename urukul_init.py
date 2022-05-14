from artiq.experiment import *
from numpy import int32, int64


class initTest(EnvExperiment):
    """init urukul"""

    def build(self):
        """
        Ensure that the necessary devices are set.
        """
        self.setattr_device("core")                                 # always needed
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")

        #self.frequency_mu = self.urukul0_ch0.frequency_to_ftw(60 * MHz)     # set frequency here
        self.frequency_mu = self.urukul0_ch0.frequency_to_ftw(110 * MHz)     # set frequency here
        self.amplitude_mu = self.urukul0_ch0.amplitude_to_asf(0.5)           # don't adjust this

    @kernel
    def run(self):
        """
        Uncomment the lines corresponding to what you want to do.
        Ensure that the correct DDS board and channel are selected.
            e.g. for board #0, channel #0, the device is urukul0_ch0
        """
        self.core.reset()
        self.urukul1_cpld.init()           # initialization for a DDS board, only needs to happen once after turning ARTIQ on
        self.core.break_realtime()
        self.urukul1_ch0.init()            # initialization for an individual DDS channel, only needs to happen once after turning ARTIQ on
        self.core.break_realtime()         # this is the time delay, put it between consecutive events
        self.urukul1_ch1.init()            # initialization for an individual DDS channel, only needs to happen once after turning ARTIQ on
        self.core.break_realtime()         # this is the time delay, put it between consecutive events
