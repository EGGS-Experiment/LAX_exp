import labrad
import numpy as np

from time import sleep
from os import environ
from artiq.experiment import *
from EGGS_labrad.config.dc_config import dc_config


class TTLTriggerVoltageSweep(EnvExperiment):
    """
    TTL Trigger Voltage Sweep
    """
    kernel_invariants = {
        'time_timeout_pmt_mu',
        'time_slack_mu',
        'time_timeout_rf_mu',
        'dc_micromotion_channels',
        'ampl_mod_vpp',
        'freq_mod_mhz',
        'dc_micromotion_voltages_v'
    }

    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge"
    ]


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # repetitions
        self.setattr_argument("repetitions",                        NumberValue(default=50000, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_us",                NumberValue(default=250, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_slack_us",                      NumberValue(default=5, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_timeout_rf_us",                 NumberValue(default=10, ndecimals=5, step=1, min=1, max=1000000))

        # modulation
        self.setattr_argument("freq_mod_mhz",                       NumberValue(default=1.4, ndecimals=5, step=1, min=0, max=1000000))
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=0.4, ndecimals=3, step=1, min=1, max=1000000))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channels",            EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list",     Scannable(
                                                                        default=RangeScan(40.0, 80.0, 41, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        # datasets
        self.set_dataset("ttl_trigger", [])
        self.setattr_dataset("ttl_trigger")

        #self.set_dataset("ttl_trigger_processed", np.zeros(self.repetitions))
        #self.setattr_dataset("ttl_trigger_processed")


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
        self.time_slack_mu =                                        self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_timeout_rf_mu =                                   self.core.seconds_to_mu(self.time_timeout_rf_us * us)

        # get voltage parameters
        self.dc_micromotion_voltages_v_list =                       np.array(list(self.dc_micromotion_voltages_v_list))
        self.dc_micromotion_channels =                              self.dc_micromotion_channeldict[self.dc_micromotion_channels]['num']

        # set up modulation
        self.fg.select_device(1)
        self.fg.frequency(self.freq_mod_mhz * 1e6)
        self.fg.toggle(1)
        self.fg.amplitude(self.ampl_mod_vpp)

        # record parameters
        self.set_dataset('xArr', self.dc_micromotion_voltages_v_list)
        self.set_dataset('repetitions', self.repetitions)
        self.set_dataset('freq_mod_mhz', self.freq_mod_mhz)
        self.set_dataset('modulation_amplitude_vpp', self.ampl_mod_vpp)
        self.set_dataset('dc_channel_num', self.dc_micromotion_channels)


    @kernel(flags='fast-math')
    def run(self):
        self.core.reset()

        # set ttl direction
        self.pmt_counter.input()
        self.rf_sync.input()
        self.core.break_realtime()


        # MAIN LOOP
        # sweep voltage
        for voltage_val in self.dc_micromotion_voltages_v_list:

            # set voltage
            self.voltage_set(self.dc_micromotion_channels, voltage_val)
            self.core.break_realtime()

            # add extra delay for ion to settle
            delay(1 * s)

            # get photon counts
            for i in range(self.repetitions):
                self.core.break_realtime()

                # wait for PMT count
                time_end_pmt_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
                time_input_pmt_mu = self.pmt_counter.timestamp_mu(time_end_pmt_mu)

                # check if event has fired
                if time_input_pmt_mu > 0:

                    # set RTIO time and add slack
                    # todo: make it now_mu
                    at_mu(time_input_pmt_mu)
                    delay_mu(self.time_slack_mu)

                    # get timestamp of RF event
                    time_end_rf_mu = self.rf_sync.gate_rising_mu(self.time_timeout_rf_mu)
                    time_input_rf_mu = self.rf_sync.timestamp_mu(time_end_rf_mu)

                    # close input gating
                    self.rf_sync.count(time_end_rf_mu)
                    self.pmt_counter.count(time_end_pmt_mu)
                    self.core.break_realtime()

                    # add data to dataset
                    self.update_dataset(voltage_val, time_input_pmt_mu, time_input_rf_mu)
                    self.core.break_realtime()


    # LABRAD FUNCTIONS
    @rpc(flags={"async"})
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        print('\tvoltage set: {}'.format(voltage_set_v))

    @rpc(flags={"async"})
    def frequency_set(self, freq_hz):
        """
        Set the RF to the desired frequency.
        """
        freq_set_hz = self.fg.frequency(freq_hz)
        print('\tfrequency set: {}'.format(freq_set_hz))

    @rpc(flags={"async"})
    def update_dataset(self, voltage_v, time_start_mu, time_stop_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('ttl_trigger', [voltage_v, self.core.mu_to_seconds(time_stop_mu - time_start_mu)])


    def analyze(self):
        # turn off modulation
        self.fg.toggle(0)

        # process data
        # ttl_trigger_tmp = np.array(self.ttl_trigger).reshape((len(self.dc_micromotion_voltages_v), self.repetitions, 2))
        # ind_arr = np.argsort(self.dc_micromotion_voltages_v)
        # ttl_trigger_tmp = ttl_trigger_tmp[ind_arr]
        #
        # self.ttl_trigger_processed = ttl_trigger_tmp
