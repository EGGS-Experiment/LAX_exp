from artiq.experiment import *
from numpy import int32, int64


class UrukulTest(EnvExperiment):
    """Urukul Test"""

    def build(self):
        #self.setattr_argument("frequency", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("amplitude", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("phase", NumberValue(ndecimals=0, step=1))
        #self.setattr_argument("channel", NumberValue(ndecimals=0, step=1))

        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch0")

        self.setattr_device("fastino0")
        self.th1 = self.urukul0_ch0.frequency_to_ftw(100*MHz)
        self.th2 = self.urukul0_ch0.amplitude_to_asf(0.4)

    @kernel
    def run(self):
        self.core.reset()
        #self.urukul0_cpld.init()
        #self.urukul0_ch0.init()
        self.urukul0_ch0.set_mu(self.th1, 0, self.th2)
        #self.core.reset()
        #self.urukul0_ch0.cfg_sw(1)
        #self.urukul0_ch0.set_att_mu(0xff)
        #self.urukul0_ch0.set_att()
        #self.urukul0_ch0.set_mu(self.th1)
        #th1 = self.urukul0_ch0.read32(0x0e)
        #th1 = self.urukul0_ch0.ftw_to_frequency(th1)
        #print(th1)
