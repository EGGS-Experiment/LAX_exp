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
        self.setattr_argument("repetitions",                        NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_us",                NumberValue(default=250, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_timeout_rf_us",                 NumberValue(default=10, ndecimals=5, step=1, min=1, max=1000000))

        # modulation
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=0.05, ndecimals=3, step=1, min=1, max=1000000))
        self.setattr_argument("freq_mod_mhz_list",                  Scannable(
                                                                        default=RangeScan(1.350, 1.400, 51, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channels",            EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltage_v",           NumberValue(default=40, ndecimals=3, step=1, min=1, max=1000000))


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.function_generator_server
        self.dc = self.cxn.dc_server


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # RF devices
        self.rf_sync =                                              self.get_device('ttl3')

        # convert time values to machine units
        self.time_timeout_pmt_mu =                                  self.core.seconds_to_mu(self.time_timeout_pmt_us * us)
        self.time_slack_mu =                                        self.core.seconds_to_mu(10 * us)
        self.time_timeout_rf_mu =                                   self.core.seconds_to_mu(self.time_timeout_rf_us * us)

        # get voltage parameters
        self.dc_micromotion_channels =                              self.dc_micromotion_channeldict[self.dc_micromotion_channels]['num']

        # reformat frequency list
        self.freq_mod_mhz_list =                                    np.array(list(self.freq_mod_mhz_list))

        # set voltage
        self.dc.voltage(self.dc_micromotion_channels, self.dc_micromotion_voltage_v)

        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("ttl_trigger",                             np.zeros([len(self.freq_mod_mhz_list) * self.repetitions, 2]))
        self.setattr_dataset("ttl_trigger")

        # record parameters
        self.set_dataset('xArr', self.freq_mod_mhz_list)
        self.set_dataset('repetitions', self.repetitions)
        self.set_dataset('modulation_amplitude_vpp', self.ampl_mod_vpp)
        self.set_dataset('dc_channel_num', self.dc_micromotion_channels)
        self.set_dataset('dc_channel_voltage', self.dc_micromotion_voltage_v)

        # set up modulation
        self.fg.select_device(1)
        self.fg.toggle(1)
        self.fg.amplitude(self.ampl_mod_vpp)


    @kernel(flags='fast-math')
    def run(self):
        self.core.reset()

        # set ttl directions
        self.pmt_counter.input()
        self.rf_sync.input()
        self.core.break_realtime()


        # MAIN LOOP
        # sweep voltage
        for freq_val_mhz in self.freq_mod_mhz_list:

            # set frequency
            self.frequency_set(freq_val_mhz * 1e6)
            self.core.break_realtime()

            # configure ttls
            self.pmt_counter._set_sensitivity(0)
            self.rf_sync._set_sensitivity(0)
            self.core.break_realtime()

            # set up counter
            counter = 0

            # get photon counts
            while counter < self.repetitions:
                #try:
                # wait for PMT count
                self.pmt_counter._set_sensitivity(1)
                delay_mu(self.time_timeout_pmt_mu)
                time_start_mu = self.pmt_counter.timestamp_mu(now_mu())

                # check if event has fired
                if time_start_mu > 0:

                    # set RTIO time
                    at_mu(time_start_mu)
                    delay_mu(self.time_slack_mu)

                    # start RF counting and stop PMT counting
                    with parallel:
                        self.pmt_counter._set_sensitivity(0)
                        self.rf_sync._set_sensitivity(1)

                    # get timestamp of RF event
                    delay_mu(self.time_timeout_rf_mu)
                    time_stop_mu = self.rf_sync.timestamp_mu(now_mu())

                    # stop RF counting
                    if time_stop_mu > 0:
                        at_mu(time_stop_mu)
                        delay_mu(self.time_slack_mu)
                        self.rf_sync._set_sensitivity(0)
                    else:
                        self.rf_sync._set_sensitivity(0)

                    # add data to dataset
                    counter += 1
                    self.update_dataset(freq_val_mhz, time_start_mu, time_stop_mu)
                    self.core.reset()

                # close gating otherwise
                else:
                    self.core.break_realtime()
                    self.pmt_counter._set_sensitivity(0)


    # LABRAD FUNCTIONS
    #@rpc(flags={"async"})
    @rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        print('\tvoltage set: {}'.format(voltage_set_v))

    #@rpc(flags={"async"})
    @rpc
    def frequency_set(self, freq_hz):
        """
        Set the RF to the desired frequency.
        """
        freq_set_hz = self.fg.frequency(freq_hz)
        print('\tfrequency set: {}'.format(freq_set_hz))

    @rpc(flags={"async"})
    def update_dataset(self, freq_hz, time_start_mu, time_stop_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self._dataset_counter += 1
        self.mutate_dataset(
            'ttl_trigger',
            self._dataset_counter,
            np.array([freq_hz, self.core.mu_to_seconds(time_stop_mu - time_start_mu)])
        )


    def analyze(self):
        # turn off modulation
        self.fg.toggle(0)

        # process data
        # ttl_trigger_tmp = np.array(self.ttl_trigger).reshape((len(self.freq_mod_mhz_list), self.repetitions, 2))
        # ind_arr = np.argsort(self.freq_mod_mhz_list)
        # ttl_trigger_tmp = ttl_trigger_tmp[ind_arr]
        #
        # self.ttl_trigger_processed = ttl_trigger_tmp
