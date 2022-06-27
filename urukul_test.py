from artiq.experiment import *
from numpy import int32, int64


class UrukulTest(EnvExperiment):
    """
    Urukul Test
    Test urukul stuff
    """

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
        self.frequency_mu = self.urukul0_ch0.frequency_to_ftw(140 * MHz)     # set frequency here
        self.amplitude_mu = self.urukul0_ch0.amplitude_to_asf(0.5)           # don't adjust this

    @kernel
    def run(self):
        """
        Uncomment the lines corresponding to what you want to do.
        Ensure that the correct DDS board and channel are selected.
            e.g. for board #0, channel #0, the device is urukul0_ch0
        """
        self.core.reset()
        #self.urukul1_ch1.set_mu(self.frequency_mu, 0, self.amplitude_mu)
        #self.urukul1_cpld.init()           # initialization for a DDS board, only needs to happen once after turning ARTIQ on
        #self.core.break_realtime()
        #self.urukul1_ch0.init()            # initialization for an individual DDS channel, only needs to happen once after turning ARTIQ on
        #self.core.break_realtime()         # this is the time delay, put it between consecutive events
        self.urukul1_ch1.cfg_sw(0)          # switches the DDS channel on/off (on=1, off=0)
        self.core.break_realtime()
        self.urukul1_ch1.set_mu(self.frequency_mu, 0, self.amplitude_mu)     # uncomment this line when you want to change frequency
        #self.urukul1_ch0.cfg_sw(1)
        self.urukul1_ch1.set_att(30 * dB)                   # set channel attenuation; has a range between 0dB and 31.5dB attenuation
        #self.urukul1_ch0.set_att(4 * dB)                     # 0dB attenuation corresponds to +2.1dBm

        #self.urukul1_ch2.cfg_sw(1)
        #self.urukul1_ch3.cfg_sw(1)
        #self.urukul1_cpld.cfg_sw(1, 1)
        #self.urukul0_cpld.cfg_sw(1, 1)
        #self.urukul1_ch3.cfg_sw(1)
        #self.urukul0_cpld.cfg_sw(2, 1)
