from artiq.experiment import *

class practice_exp(EnvExperiment):
    """
    testing
    """
    # def __init__(self):
    #     print(self._HasEnvironment__device_mgr)
    #
    # def build(self):
    #     self.setattr_device('core')
    #     self.setattr_device('core_dma')
    #     self.setattr_device('scheduler')
    #     self.setattr_device('ttl4')
    #     self.setattr_device('ttl5')
    #
    # @kernel
    # def _record(self):
    #     with self.core_dma.record('ps'):
    #         for i in range(400):
    #             with parallel:
    #                 self.ttl4.pulse(1*ms)
    #                 #self.ttl5.pulse(1*ms)
    #             delay(1.0*ms)
    #
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     self._record()
    #     handle = self.core_dma.get_handle('ps')
    #     self.core.reset()
    #     for i in range(100):
    #         self.core_dma.playback_handle(handle)
    def build(self):
        self.setattr_device('core')
        #self.setattr_device('core_dma')
        #self.setattr_device('scheduler')
        self.setattr_device('ttl4')
        self.setattr_device('urukul0_ch0')
        self.setattr_device('urukul0_cpld')
        self.setattr_device('urukul1_ch0')
        self.setattr_device('urukul1_cpld')
        self.setattr_device('urukul2_ch0')
        self.setattr_device('urukul2_cpld')
        self.setattr_device('zotino0')

    @kernel
    def run(self):
        self.core.reset()
        self.zotino0.init()
        self.core.break_realtime()
        self.zotino0.write_dac_mu(0, 20000)
        self.zotino0.load()
        self.core.break_realtime()
        self.zotino0.load()
        self.core.reset()
        self.urukul1_cpld.init()
        self.urukul1_ch0.init()
        self.urukul1_ch0.set_att(0*dB)
        print(self.urukul1_cpld.get_att_mu())
        self.core.reset()
        self.urukul1_ch0.cfg_sw(False)
        self.urukul1_ch0.set_amplitude(1.0)
        self.urukul1_ch0.set_frequency(20000)
        self.core.reset()
        for i in range(100):
            self.ttl4.pulse(1*ms)
            delay(1*ms)
            #self.core.close()
            #print(self.scheduler.get_status())
        # self.core.reset()
        # self.urukul2_cpld.init()
        # self.urukul2_ch0.init()
        # handle = self.core_dma.get_handle('ps')
        # self.core.reset()
        # for i in range(200):
        #     self.core_dma.playback_handle(handle)
        #     self.core.reset()

