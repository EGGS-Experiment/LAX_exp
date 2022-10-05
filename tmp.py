import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG


class Testing(EnvExperiment):
    """
    tmp exp
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
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

    def prepare(self):
        self.asf = self.urukul0_ch0.amplitude_to_asf(0.5)
        self.asf2 = self.urukul0_ch0.amplitude_to_asf(0.2)
        self.ftw = self.urukul0_ch0.frequency_to_ftw(50 * MHz)
        self.ftw2 = self.urukul0_ch0.frequency_to_ftw(50 * MHz)
        # self.amplitude = 0
        # self.frequency = 0
        # self.phase = 0
        # self.regval = np.int64(0)
        # self.regval2 = 0
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.urukul0_cpld.cfg_switches(0b1111)
        self.urukul0_cpld.set_profile(1)
        self.urukul0_cpld.cfg_switches(0b1111)
        #self.urukul0_ch3.set_mu(self.ftw2, asf=self.asf2, profile=1)
        # self.core.break_realtime()
        # self.urukul0_ch3.set_mu(self.ftw, asf=self.asf, profile=0)
        #self.core.break_realtime()




    def analyze(self):
        pass
        # ftw = np.int32(self.regval & 0xffffffff)
        # ampl = np.int32((self.regval >> 48) & 0xffff)
        # print('ftw:', ftw)
        # print('ampl:', ampl)
