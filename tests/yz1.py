import labrad
import numpy as np

from os import environ
from artiq.experiment import *
from EGGS_labrad.config.dc_config import dc_config


class TTLTriggerFrequencySweep(EnvExperiment):
    """
    TTL Trigger Frequency Sweep
    """
    kernel_invariants = {
        'time_timeout_pmt_mu',
        'time_slack_mu',
        'time_timeout_rf_mu',
        'dc_micromotion_channels',
        'ampl_mod_vpp',
        'freq_mod_mhz_list',
        'dc_micromotion_voltage_v'
    }

    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge"
    ]


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # repetitions
        self.setattr_argument("repetitions",                        NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_us",                NumberValue(default=25000, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_slack_us",                      NumberValue(default=5, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_timeout_rf_us",                 NumberValue(default=10, ndecimals=5, step=1, min=1, max=1000000))

        # datasets
        self.set_dataset("ttl_trigger", np.zeros([self.repetitions, 2]))
        self.setattr_dataset("ttl_trigger")

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.function_generator_server


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # RF devices
        self.rf_sync =                                              self.get_device('ttl3')

        # convert time values to machine units
        self.time_timeout_pmt_mu =                                  self.core.seconds_to_mu(self.time_timeout_pmt_us * us)
        self.time_slack_mu =                                        self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_timeout_rf_mu =                                   self.core.seconds_to_mu(self.time_timeout_rf_us * us)

        # general
        self.loop_iter = list(range(self.repetitions))

        # set up modulation
        self.fg.select_device(1)
        self.fg.toggle(1)
        self.fg.frequency(1400000)

    @kernel(flags='fast-math')
    def run(self):
        self.core.reset()

        # set ttl direction
        #time_start_mu = 0
        #time_stop_mu = 0
        self.core.break_realtime()


        # MAIN LOOP
        for i in self.loop_iter:
            self.core.reset()

            # start PMT counting
            self.rf_sync._set_sensitivity(0)
            self.pmt_counter._set_sensitivity(1)
            #delay_mu(self.time_timeout_pmt_mu)
            delay_mu(100000)
            time_start_mu = self.pmt_counter.timestamp_mu(now_mu())

            if time_start_mu > 0:
                at_mu(time_start_mu)
                # DELAY IS NEEDED
                delay_mu(100)

                # stop PMT counting, start RF counting
                self.pmt_counter._set_sensitivity(0)
                self.rf_sync._set_sensitivity(1)

                # stop RF counting
                time_stop_mu = self.rf_sync.timestamp_mu(now_mu())

                at_mu(time_stop_mu)
                delay_mu(100)
                self.rf_sync._set_sensitivity(0)

                self.core.break_realtime()
                #self.update_dataset(i, time_start_mu, time_stop_mu)
                self.core.break_realtime()
            self.core.reset()


    #
    #
    # @kernel(flags='fast-math')
    # def DMArecord(self):
    #     with self.core_dma.record('yzde'):
    #
    #         # RESET
    #             # set cooling waveform
    #         with parallel:
    #             self.dds_board.set_profile(0)
    #             delay_mu(self.time_profileswitch_delay_mu)
    #
    #             # qubit repump
    #         self.dds_board.cfg_switches(0b1100)
    #         delay_mu(self.time_qubit_repump_ini
    @rpc(flags={"async"})
    def update_dataset(self, i, time_start_mu, time_stop_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.mutate_dataset('ttl_trigger', i, [100, 100])


    def analyze(self):
        print(self.ttl_trigger)
        # turn off modulation
        #self.fg.toggle(0)

        # process data
        # ttl_trigger_tmp = np.array(self.ttl_trigger).reshape((len(self.freq_mod_mhz_list), self.repetitions, 2))
        # ind_arr = np.argsort(self.freq_mod_mhz_list)
        # ttl_trigger_tmp = ttl_trigger_tmp[ind_arr]
        #
        # self.ttl_trigger_processed = ttl_trigger_tmp
