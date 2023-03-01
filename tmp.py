import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testarg12(EnvExperiment):
    """
    testarg34
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        self.setattr_device("ttl10")    # rf blank
        self.setattr_device("ttl11")    # oscilloscope trigger
        self.setattr_device("ttl12")    # rf switch

        # self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        # self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
        #                                                                     default=CenterScan(85, 5, 0.2, randomize=True),
        #                                                                     global_min=30, global_max=400, global_step=1,
        #                                                                     unit="MHz", scale=1, ndecimals=5
        #                                                                 ))

        # self.setattr_device("core_dma")
        # self.setattr_device('urukul1_ch0')

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    def prepare(self):
        # arg vals
        self.num_counts =           1
        self.pwm_duty_cycle_pct =   50
        self.pwm_freq_hz =          100000
        self.time_holdoff_ns =      10

        # tmp remove
        self.time_pulseoff_us =     20
        self.time_pulseoff_mu =     self.core.seconds_to_mu(self.time_pulseoff_us * us)

        # exp vals
        self._iter_loop =           np.arange(self.num_counts)
        self._handle_dma =          "TMP_DMA"

        # conv vals
        self.pwm_delay_on_mu =      self.core.seconds_to_mu((self.pwm_duty_cycle_pct / 100) / self.pwm_freq_hz)
        self.pwm_delay_off_mu =     self.core.seconds_to_mu((1 - self.pwm_duty_cycle_pct / 100) / self.pwm_freq_hz)
        self.time_holdoff_mu =      self.core.seconds_to_mu(self.time_holdoff_ns * ns)

        # alias devices
        self.ttl_rf =               self.get_device('ttl10')
        self.ttl_os =               self.get_device('ttl11')
        self.ttl_sw =               self.get_device('ttl12')


    @kernel
    def run(self):
        # todo: record dma

        # reset
        self.core.reset()
        # with parallel:
        self.ttl_rf.off()
        self.ttl_os.off()
        self.ttl_sw.off()
        delay_mu(5000)

        # trigger scope
        self.ttl_os.on()
        delay_mu(self.time_holdoff_mu)

        # # blank PWM
        # for i in self._iter_loop:
        #     self.ttl_sw.on()
        #     delay_mu(self.pwm_delay_on_mu)
        #     self.ttl_sw.off()
        #     delay_mu(self.pwm_delay_off_mu)

        # one-shot pulse
        self.ttl_sw.on()
        delay_mu(self.time_pulseoff_mu)
        self.ttl_sw.off()
        delay_mu(self.time_pulseoff_mu)

        # unset
        self.ttl_rf.off()
        self.ttl_os.off()
        self.ttl_sw.off()

    @kernel
    def recordDMA(self):
        pass

    def analyze(self):
        print('test done')
