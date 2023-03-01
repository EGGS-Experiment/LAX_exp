import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class trapPWM(EnvExperiment):
    """
    testarg34
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        self.setattr_device("ttl10")    # rf blank
        self.setattr_device("ttl11")    # oscilloscope trigger
        self.setattr_device("ttl12")    # rf switch - main

        # pwm
        self.pwm_freq_hz =                                      10

        # timing
        self.time_pwm_off_us =                                  3
        self.time_holdoff_ns =                                  10
        self.time_run_s =                                       5


    def prepare(self):
        # exp vals
        self.time_rfswitch_delay_mu =                           self.core.seconds_to_mu(2 * us)
        self._iter_loop =                                       np.arange(int(self.time_run_s * (self.pwm_freq_hz * 1)))
        self._handle_dma =                                      "TMP_DMA"

        # conv vals
        self.time_pwm_off_mu =                                  self.core.seconds_to_mu(self.time_pwm_off_us * us)
        self.time_pwm_on_mu =                                   self.core.seconds_to_mu(1 / (self.pwm_freq_hz * 1) - self.time_pwm_off_us * us)
        self.time_holdoff_mu =                                  self.core.seconds_to_mu(self.time_holdoff_ns * ns)

        # alias devices
        self.ttl_rf =                                           self.get_device('ttl10')
        self.ttl_os =                                           self.get_device('ttl11')
        self.ttl_sw_main =                                      self.get_device('ttl12')
        self.ttl_sw_cancel =                                    self.get_device('ttl13')


    @kernel
    def run(self):
        # reset
        self.core.reset()
        with parallel:
            self.ttl_rf.off()
            self.ttl_os.off()
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()
        delay_mu(5000)

        # record dma and get handle
        self.DMArecord()
        _handle_tmp = self.core_dma.get_handle(self._handle_dma)
        self.core.break_realtime()

        # MAIN LOOP
        for i in self._iter_loop:
            self.core_dma.playback_handle(_handle_tmp)

        # finish and unset TTLs
        self.core.break_realtime()
        with parallel:
            self.ttl_rf.off()
            self.ttl_os.off()
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # PWM sequence
        with self.core_dma.record(self._handle_dma):

            # rf off
            with parallel:
                self.ttl_os.on()
                self.ttl_sw_main.on()
            delay_mu(self.time_pwm_off_mu)

            # rf on
            with parallel:
                self.ttl_os.off()
                self.ttl_sw_main.off()
                delay_mu(self.time_pwm_on_mu)


    def analyze(self):
        print('test done')
