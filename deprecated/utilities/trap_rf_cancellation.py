import numpy as np
from artiq.experiment import *
_DMA_HANDLE_SEQUENCE = "TRAP_PWM_SEQUENCE"


class TrapRFCancellation(EnvExperiment):
    """
    Trap RF Cancellation

    Pulses the trap RF in response to a trigger to minimize the electric field during loading.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("time_run_s",                     NumberValue(default=10, ndecimals=5, step=0.1, min=0, max=1000000))

        # trap cancellation parameters
        self.setattr_argument("time_rf_cancel_us",              NumberValue(default=0.5, ndecimals=5, step=0.1, min=0, max=1000000))
        self.setattr_argument("time_rf_restore_us",             NumberValue(default=0.5, ndecimals=5, step=0.1, min=0, max=1000000))

        # trigger timing
        self.setattr_argument("freq_trig_hz",                   NumberValue(default=10, ndecimals=3, step=0.1, min=0.001, max=10000))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # alias devices
        self.counter_trigger =                                  self.get_device('ttl4')
        self.counter_rf =                                       self.get_device('ttl3')
        self.ttl_sw_main =                                      self.get_device('ttl12')
        self.ttl_sw_cancel =                                    self.get_device('ttl13')

        # experiment runs
        self._iter_loop =                                       np.arange(int(self.time_run_s * (self.freq_trig_hz * 1)))

        # trap cancellation values
        self.time_rf_cancel_mu =                                self.core.seconds_to_mu(self.time_rf_cancel_us * us)
        self.time_rf_restore_mu =                               self.core.seconds_to_mu(self.time_rf_restore_us * us)

        # trigger values
        self.time_trigger_gating_mu =                           self.core.seconds_to_mu(1 / (self.freq_trig_hz * 1))
        self.time_trigger_delay_mu =                            self.core.seconds_to_mu(50.88 * us)


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        _handle_tmp = self.core_dma.get_handle(_DMA_HANDLE_SEQUENCE)
        self.core.break_realtime()

        # tmp remove
        tmphandle = self.core_dma.get_handle('tmpseq')


        # MAIN LOOP
        for i in self._iter_loop:
            self.core.break_realtime()

            # wait for trigger input
            time_trigger_mu = self.counter_trigger.timestamp_mu(self.counter_trigger.gate_rising_mu(self.time_trigger_gating_mu))

            # respond to trigger input
            if time_trigger_mu > 0:

                # set rtio clock to input trigger time
                at_mu(time_trigger_mu + 1000)
                self.core.break_realtime()

                # wait until next RF pulse
                time_rf_mu = self.counter_rf.timestamp_mu(self.counter_rf.gate_rising_mu(1000))

                if time_rf_mu > 0:
                    #print(self.core.get_rtio_counter_mu() - now_mu())
                    #self.core.break_realtime()
                    # delay_mu(100000)
                    # print(now_mu() - time_rf_mu)
                    # at_mu(time_rf_mu + 1000)
                    #self.core.break_realtime()

                    at_mu(time_rf_mu + 12000)
                    #with parallel:
                    self.counter_rf._set_sensitivity(0)
                    #at_mu(self.core.get_rtio_counter_mu() + 500000)

                    # run pulse sequence
                    #self.core.break_realtime()
                    self.core_dma.playback_handle(_handle_tmp)
                    self.core.break_realtime()

            # wait for next pulse
            else:
                self.core.break_realtime()

            # tmp remove
            self.core.reset()


        # finish and unset TTLs
        self.core.break_realtime()
        with parallel:
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # trap sequence
        with self.core_dma.record(_DMA_HANDLE_SEQUENCE):

            # match trigger delay
            delay_mu(self.time_trigger_delay_mu)

            # rf off
            with parallel:
                self.ttl_sw_main.on()
                self.ttl_sw_cancel.off()

            # actively cancel residual rf
            delay_mu(self.time_rf_cancel_mu)

            # restore residual rf
            self.ttl_sw_cancel.on()
            delay_mu(self.time_rf_restore_mu)

            # rf on
            with parallel:
                self.ttl_sw_main.off()
                self.ttl_sw_cancel.off()

        # tmp remove
        # trap sequence
        with self.core_dma.record('tmpseq'):
            with parallel:
                self.counter_trigger._set_sensitivity(0)
                self.counter_rf._set_sensitivity(1)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # prepare pulse sequence TTL triggers
        with parallel:
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()


    def analyze(self):
        print('test done')
