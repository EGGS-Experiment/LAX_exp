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
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")
        self.setattr_device("ttl10")

        self.time0 = self.core.seconds_to_mu(1000 * ms)
        self.time1 = self.core.seconds_to_mu(100 * us)

    #@kernel
    def run(self):
        self.set_dataset('dds_tickle_channel', 3, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_qubit_pct', 14.0, broadcast=True, persist=True)
        #
        # self.set_dataset('freq_pump_cooling_mhz', 110.0, broadcast=True, persist=True)
        # self.set_dataset('freq_pump_readout_mhz', 110.0, broadcast=True, persist=True)
        # self.set_dataset('freq_repump_cooling_mhz', 110.0, broadcast=True, persist=True)
        # self.set_dataset('freq_repump_qubit_mhz', 110.0, broadcast=True, persist=True)
       # self.set_dataset('ttl_channel_function_generator', 9, broadcast=True, persist=True)
    #     self.core.reset()
    #     self.ttl8.off()
    #     self.ttl9.off()
    #     self.ttl10.off()
    #     self.core.reset()
    #
    #     delay_mu(self.time0)
    #
    #     self.ttl8.on()
    #     delay_mu(self.time1)
    #
    #     with parallel:
    #         self.ttl9.on()
    #         self.ttl10.on()
    #
    #     delay_mu(self.time1)
    #
    #     with parallel:
    #         self.ttl8.off()
    #         self.ttl9.off()
    #         self.ttl10.off()        #pass
    # #
    # # def run(self):
    #     pass
