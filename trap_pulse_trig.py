import labrad
import numpy as np
from os import environ
from artiq.experiment import *
# todo: examine pulse delay on scope and account


class TrapPulseTrigger(EnvExperiment):
    """
    Trap Pulse Trigger

    todo add description
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        self.setattr_device("ttl10")                            # rf blank
        self.setattr_device("ttl11")                            # oscilloscope trigger
        self.setattr_device("ttl12")                            # rf switch - main
        self.setattr_device("ttl13")                            # rf switch - cancellation
        self.setattr_device("ttl4")                             # 243nm trigger - input

        # general timing
        self.time_run_s =                                       600

        # pwm timing
        self.time_pwm_cancel_us =                               0.5
        self.time_pwm_restore_us =                              0.5

        # 243 trigger timing
        self.freq_trig_243_hz =                                 10
        self.time_delay_243_us =                                50.88
        self.time_243_holdoff_ns =                              10


    def prepare(self):
        # alias devices
        self.ttl_rf =                                           self.get_device('ttl10')
        self.ttl_os =                                           self.get_device('ttl11')
        self.ttl_sw_main =                                      self.get_device('ttl12')
        self.ttl_sw_cancel =                                    self.get_device('ttl13')
        self.counter_trig_243 =                                 self.get_device('ttl4')

        # exp vals
        self._iter_loop =                                       np.arange(int(self.time_run_s * (self.freq_trig_243_hz * 1)))
        self._handle_dma =                                      "TRAP_PWM_SEQUENCE"

        # conv pwm vals
        self.time_pwm_cancel_mu =                               self.core.seconds_to_mu(self.time_pwm_cancel_us * us)
        self.time_pwm_restore_mu =                              self.core.seconds_to_mu(self.time_pwm_restore_us * us)

        # conv 243 vals
        self.time_trigger_243_gating_mu =                       self.core.seconds_to_mu(1 / (self.freq_trig_243_hz * 1))
        self.time_delay_243_mu =                                self.core.seconds_to_mu(self.time_delay_243_us * us)
        self.time_243_holdoff_mu =                              self.core.seconds_to_mu(self.time_243_holdoff_ns * ns)


    @kernel
    def run(self):
        # reset
        self.core.reset()
        with parallel:
            self.ttl_rf.off()
            self.ttl_os.off()
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()
            delay_mu(1000)

        # record dma and get handle
        self.DMArecord()
        _handle_tmp = self.core_dma.get_handle(self._handle_dma)
        self.core.break_realtime()


        # MAIN LOOP
        for i in self._iter_loop:
            self.core.break_realtime()

            # wait for trigger input
            time_trigger_mu = self.counter_trig_243.timestamp_mu(self.counter_trig_243.gate_rising_mu(self.time_trigger_243_gating_mu))

            # respond to trigger input
            if time_trigger_mu > 0:

                # set rtio clock to input trigger time
                at_mu(time_trigger_mu)
                self.core.break_realtime()

                # run pulse sequence
                self.core_dma.playback_handle(_handle_tmp)

            # wait for next pulse
            else:
                self.core.break_realtime()


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

            # match 243 trigger time
            delay_mu(self.time_delay_243_mu)

            # rf off
            with parallel:
                self.ttl_os.on()
                self.ttl_sw_main.on()
                self.ttl_sw_cancel.off()

            # actively cancel residual rf
            delay_mu(self.time_pwm_cancel_mu)

            # restore residual rf
            self.ttl_sw_cancel.on()
            delay_mu(self.time_pwm_restore_mu)

            # rf on
            with parallel:
                self.ttl_os.off()
                self.ttl_sw_main.off()
                self.ttl_sw_cancel.off()


    def analyze(self):
        print('test done')
