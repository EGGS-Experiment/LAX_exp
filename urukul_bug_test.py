from artiq.experiment import *
from artiq.coredevice.ad9910 import _AD9910_REG_FTW

class urukul_exp(EnvExperiment):
    """
    testing urukul
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_device('ttl4')
        self.setattr_device('urukul0_ch0')
        self.setattr_device('urukul0_ch1')
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_ch0')
        self.setattr_device('urukul1_cpld')
        print('scde:', type(self.urukul1_ch0))
        print('scde:', self.urukul1_ch0.__class__)
        print('scde:', self.urukul1_ch0.__class__.__name__)
        #self.setattr_device('urukul2_ch0')
        #self.setattr_device('urukul2_cpld')
        #self.setattr_device('zotino0')
        self.setattr_device('fastino0')
        self.freq1=200*MHz
        self.yz1 = 0

    @kernel
    def run(self):
        self.core.reset()
        # self.urukul0_ch0.write32(0x07, 0xffffffff)
        # self.core.break_realtime()
        # # #self.urukul0_ch0.set_ftw(0xffffffff)
        # rx1 = self.urukul0_ch0.get_mu()
        # self.core.reset()
        # print(rx1)
        # self.core.reset()
        # self.yz1 = self.urukul0_ch0.read32(0x07)
        # self.core.reset()
        # print(self.yz1)
        # self.urukul1_cpld.init()
        # self.core.reset()
        # self.urukul1_ch0.init()
        #self.urukul0_cpld.set_att(0, 0*dB)
        #self.urukul1_cpld.init()
        self.core.reset()
        self.urukul0_cpld.init()
        self.core.reset()
        self.urukul1_cpld.init()
        self.urukul0_ch0.cfg_sw(True)
        self.core.reset()
        self.fastino0.init()




