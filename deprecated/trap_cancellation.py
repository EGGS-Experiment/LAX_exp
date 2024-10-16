import labrad
import numpy as np
from os import environ
from time import sleep
from numpy import int32, int64

from artiq.experiment import *
from artiq.coredevice.rtio import *

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
        self.setattr_argument("time_run_s",                         NumberValue(default=10, precision=5, step=0.1, min=0, max=1000000))

        # trap cancellation parameters
        self.setattr_argument("time_rf_cancel_us",                  NumberValue(default=0.5, precision=5, step=0.1, min=0, max=1000000))
        self.setattr_argument("time_rf_restore_us",                 NumberValue(default=0.5, precision=5, step=0.1, min=0, max=1000000))

        # trigger timing
        self.setattr_argument("freq_trig_hz",                       NumberValue(default=10, precision=3, step=0.1, min=0.001, max=10000))
        self.setattr_argument("freq_rf_mhz",                        NumberValue(default=19.0435, precision=7, step=0.001, min=0.000001, max=10000000))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # alias devices
        self.counter_trigger =                                      self.get_device('ttl4')
        self.counter_rf =                                           self.get_device('ttl7')
        self.ttl_sw_main =                                          self.get_device('ttl12')
        self.ttl_sw_cancel =                                        self.get_device('ttl13')
        self.ttl_trig_cancel =                                      self.get_device('ttl8')

        # experiment runs
        self._iter_loop =                                           np.arange(int(self.time_run_s * (self.freq_trig_hz * 1)))

        # trigger values
        self.time_trigger_gating_mu =                               self.core.seconds_to_mu(1 / (self.freq_trig_hz * 1))
        self.time_trigger_delay_mu =                                self.core.seconds_to_mu(161.38 * us)

        # rf values (for synchronization)
        self.time_rf_period_s =                                     1 / (self.freq_rf_mhz * MHz)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(2 * self.time_rf_period_s)

        # try to sync cancellation time to rf period
        self.num_periods_cancel =                                   round((self.time_rf_cancel_us * us) / self.time_rf_period_s)
        self.num_periods_restore =                                  round((self.time_rf_restore_us * us) / self.time_rf_period_s)
        self.time_rf_cancel_mu =                                    self.core.seconds_to_mu(self.num_periods_cancel * self.time_rf_period_s)
        self.time_rf_restore_mu =                                   self.core.seconds_to_mu(self.num_periods_restore * self.time_rf_period_s)

        # tmp remove
        self._iter_dataset =                                        0
        self.set_dataset("rf_cancellation",                         np.zeros(len(self._iter_loop)))
        self.setattr_dataset("rf_cancellation")


        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg =                                                   self.cxn.function_generator_server

        # get list of function generators
        fg_dev_list = self.fg.list_devices()
        fg_dev_dict = dict(tuple(fg_dev_list))

        # select correct function generator
        dev_exists = False
        for dev_num, dev_desc in fg_dev_dict.items():
            if 'DG2P' in dev_desc:
                dev_exists = True
                self.fg.select_device(dev_num)

        # raise error if function generator doesn't exist
        if not dev_exists:
            raise Exception("Error: modulation function generator not detected.")


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
        _handle_seq = self.core_dma.get_handle(_DMA_HANDLE_SEQUENCE)
        self.core.break_realtime()


        # MAIN LOOP
        for i in self._iter_loop:

            # wait for trigger input
            self.counter_trigger._set_sensitivity(1)
            time_trigger_mu = self.counter_trigger.timestamp_mu(now_mu() + self.time_trigger_gating_mu)

            # respond to trigger input
            if time_trigger_mu > 0:

                # set rtio hardware time to input trigger time
                at_mu(time_trigger_mu + 7000)
                # wait until next RF pulse for synchronization
                self.counter_rf._set_sensitivity(1)
                time_rf_mu = self.counter_rf.timestamp_mu(now_mu() + self.time_rf_gating_mu)

                # start RF sequence
                if time_rf_mu > 0:
                    # set rtio hardware time to RF pulse time
                    at_mu(time_rf_mu + 8000)
                    # run pulse sequence
                    self.core_dma.playback_handle(_handle_seq)

                    # # save timing data
                    self.core.break_realtime()
                    self.update_dataset(time_rf_mu - time_trigger_mu)
                    self.core.break_realtime()

            # reset input buffers
            self.core.reset()


        # finish and unset TTLs
        self.core.break_realtime()
        with parallel:
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()
            self.ttl_trig_cancel.off()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # trap sequence
        with self.core_dma.record(_DMA_HANDLE_SEQUENCE):

            # stop rf counting and match trigger delay
            with parallel:
                self.ttl_trig_cancel.on()
                self.counter_rf._set_sensitivity(0)
                self.counter_trigger._set_sensitivity(0)
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
                self.ttl_trig_cancel.off()


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set up function generator
        self.fg.gpib_write(':OUTP:IMP 50')
        self.fg.toggle(0)
        self.fg.burst(True)
        self.fg.burst_mode('GAT')
        # todo: set amplitude
        # self.fg.amplitude(self.ampl_mod_vpp)
        # todo: couple CH1 and CH2 triggering
        self.fg.gpib_write(':COUP1:TRI 1')
        sleep(0.2)
        # todo: phase lock and align CH1 and CH2
        self.fg.gpib_write(':COUP1:PHAS:MODE OFFS')
        sleep(0.2)
        self.fg.gpib_write(':COUP1:PHAS:STAT 1')
        sleep(0.2)
        self.fg.gpib_write(':COUP2:PHAS:DEV 90')
        sleep(0.2)
        # todo: set output phase
        # self.fg.phase(0)
        self.fg.toggle(1)

        # prepare pulse sequence TTL triggers
        with parallel:
            self.ttl_sw_main.off()
            self.ttl_sw_cancel.off()
            self.ttl_trig_cancel.off()


    @rpc(flags={"async"})
    def update_dataset(self, val):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # save data to datasets and update iterator
        self.mutate_dataset('rf_cancellation', self._iter_dataset, val)
        self._iter_dataset += 1


    def analyze(self):
        print('test done')
        print('\ttimings: {:.4f} +/- {:.4f}'.format(np.mean(self.rf_cancellation), np.std(self.rf_cancellation)))
        print('\trange: {:.4f}'.format(np.max(self.rf_cancellation) - np.min(self.rf_cancellation)))
        print('\tmax: {:.4f}, min: {:.4f}'.format(np.max(self.rf_cancellation), np.min(self.rf_cancellation)))
