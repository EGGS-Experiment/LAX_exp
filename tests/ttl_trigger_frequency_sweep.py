import labrad
import numpy as np

from time import sleep
from os import environ
from artiq.experiment import *
from EGGS_labrad.config.dc_config import dc_config


class TTLTriggerFrequencySweep(EnvExperiment):
    """
    TTL Trigger Frequency Sweep
    """
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge"
    ]


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # repetitions
        self.setattr_argument("repetitions",                        NumberValue(default=20000, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_us",                NumberValue(default=1000, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_slack_us",                      NumberValue(default=5, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_timeout_rf_us",                 NumberValue(default=10, ndecimals=5, step=1, min=1, max=1000000))

        # modulation
        self.setattr_argument("freq_mod_mhz_list",                  Scannable(
                                                                        default=RangeScan(1.470, 1.500, 31),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=2.0, ndecimals=3, step=1, min=1, max=1000000))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channels",            EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='A-Ramp 1'))
        self.setattr_argument("dc_micromotion_voltage_v",           NumberValue(default=40, ndecimals=3, step=1, min=1, max=1000000))

        # datasets
        self.set_dataset("ttl_trigger", [])
        self.setattr_dataset("ttl_trigger")

        self.set_dataset("ttl_trigger_processed", np.zeros(self.repetitions))
        self.setattr_dataset("ttl_trigger_processed")


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.fg_server
        self.dc = self.cxn.dc_server


    def prepare(self):
        # PMT devices
        self.pmt_counter =                      self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # RF devices
        self.rf_sync =                          self.get_device('ttl3')

        # convert time values to machine units
        self.time_timeout_pmt_mu =              self.core.seconds_to_mu(self.time_timeout_pmt_us * us)
        self.time_slack_mu =                    self.core.seconds_to_mu(self.time_slack_us * us)
        self.time_timeout_rf_mu =               self.core.seconds_to_mu(self.time_timeout_rf_us * us)

        # get voltage parameters
        self.dc_micromotion_channels =          self.dc_micromotion_channeldict[self.dc_micromotion_channels]['num']

        # set voltage
        self.dc.voltage(self.dc_micromotion_channels, self.dc_micromotion_voltage_v)

        # set up modulation
        self.fg.select_device()
        self.fg.gpib_write('OUTP ON')
        self.fg.gpib_write('VOLT {}'.format(self.ampl_mod_vpp))

        # tmp remove: record parameters
        self.set_dataset('xArr', self.frequency_list_mhz)
        self.set_dataset('dc_channel_num', self.dc_micromotion_channels)
        self.set_dataset('dc_channel_voltage', self.dc_micromotion_voltage_v)


    @kernel
    def run(self):
        self.core.reset()


        # MAIN LOOP
        # sweep voltage
        for freq_val_mhz in self.freq_mod_mhz_list:

            # set frequency
            self.frequency_set(freq_val_mhz * 1e6)
            self.core.break_realtime()

            # get photon counts
            for i in range(self.repetitions):
                self.core.break_realtime()

                # wait for PMT count
                time_end_pmt_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
                time_input_pmt_mu = self.pmt_counter.timestamp_mu(time_end_pmt_mu)

                # check if event has fired
                if time_input_pmt_mu > 0:

                    # set RTIO time and add slack
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
                    self.update_dataset(freq_val_mhz, time_end_pmt_mu, time_input_rf_mu)
                    self.core.break_realtime()

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
        self.fg.gpib_write('FREQ {}'.format(freq_hz))
        sleep(1.0)
        freq_set_hz = self.fg.gpib_query('FREQ?')
        print('frequency set: {}'.format(freq_set_hz))

    @rpc(flags={"async"})
    def update_dataset(self, freq_hz, time_start_mu, time_stop_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('ttl_trigger', [freq_hz, time_start_mu, time_stop_mu])


    def analyze(self):
        ttl_trigger_tmp = np.array(self.ttl_trigger).reshape((len(self.dc_micromotion_voltages_v, self.repetitions, 2)))
        ind_arr = np.argsort(self.freq_mod_mhz_list)
        ttl_trigger_tmp = ttl_trigger_tmp[ind_arr]

        self.ttl_trigger_processed = ttl_trigger_tmp
