import labrad
import numpy as np
from time import sleep

from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from EGGS_labrad.config.dc_config import dc_config


class ParametricSweep(EnvExperiment):
    """
    Parametric Sweep
    """
    kernel_invariants = {
        'time_pmt_gating_mu',
        'dc_micromotion_channel',
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

        # num_counts
        self.setattr_argument("num_counts",                         NumberValue(default=100000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=0.2, ndecimals=3, step=0.01, min=0, max=1000000))
        self.setattr_argument("freq_mod_mhz_list",                  Scannable(
                                                                        default=CenterScan(1.385, 0.01, 0.0002, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=0.001,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel",             EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltage_v",           NumberValue(default=50.0, ndecimals=3, step=1, min=0, max=1000000))


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))
        self.time_pmt_gating_mu =                                   self.core.seconds_to_mu(100 * us)

        # get voltage parameters
        self.set_dataset('dc_channel_name', self.dc_micromotion_channel)
        self.dc_micromotion_channel =                               self.dc_micromotion_channeldict[self.dc_micromotion_channel]['num']

        # reformat frequency list
        self.freq_mod_mhz_list =                                    np.array(list(self.freq_mod_mhz_list))

        # modulation control and synchronization
        self.mod_toggle =                                           self.get_device("ttl8")
        self.mod_clock =                                            self.get_device("urukul0_ch3")
        self.mod_clock_freq_ftw =                                   self.mod_clock.frequency_to_ftw(10. * MHz)
        self.mod_clock_ampl_pct =                                   self.mod_clock.amplitude_to_asf(0.5)
        self.mod_clock_att_db =                                     4 * dB
        # self.time_mod_delay_mu =                                    self.core.seconds_to_mu(300 * ns)
        self.time_mod_delay_mu =                                    self.core.seconds_to_mu(300 * ns)

        # RF synchronization
        self.rf_clock =                                             self.get_device('ttl7')
        self.time_rf_holdoff_mu =                                   self.core.seconds_to_mu(10000 * ns)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(100 * ns)


        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([len(self.freq_mod_mhz_list), self.num_counts]))
        self.setattr_dataset("results")

        # record parameters
        self.set_dataset('xArr',                                    self.freq_mod_mhz_list)
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_amplitude_vpp',                self.ampl_mod_vpp)
        self.set_dataset('dc_channel_n',                            self.dc_micromotion_channel)
        self.set_dataset('dc_channel_voltage',                      self.dc_micromotion_voltage_v)

        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg =                                                   self.cxn.function_generator_server
        self.dc =                                                   self.cxn.dc_server

        # set up function generator
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

        # set up function generator
        self.fg.gpib_write(':OUTP:IMP 50')
        self.fg.toggle(0)
        self.fg.amplitude(self.ampl_mod_vpp)
        self.fg.burst(True)
        self.fg.burst_mode('GAT')
        self.fg.toggle(1)

        # set voltage
        #self.dc.voltage(self.dc_micromotion_channel, self.dc_micromotion_voltage_v)


    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device
        self.core.reset()

        # set ttl directions
        with parallel:
            self.pmt_counter.input()
            self.rf_clock.input()
            self.mod_toggle.output()
            self.mod_toggle.off()
        self.core.break_realtime()

        # configure rf mod clock
        self.mod_clock.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.mod_clock.set_att(self.mod_clock_att_db)
        self.mod_clock.set_mu(self.mod_clock_freq_ftw, asf=self.mod_clock_ampl_pct)
        self.mod_clock.cfg_sw(True)
        self.core.break_realtime()

        # enable external clocking
        self.fg_write(':ROSC:SOUR EXT')
        self.core.break_realtime()

        # set voltage
        self.voltage_set(self.dc_micromotion_channel, self.dc_micromotion_voltage_v)
        self.core.break_realtime()


        # MAIN LOOP
        # sweep modulation frequency
        for freq_val_mhz in self.freq_mod_mhz_list:

            self.core.break_realtime()

            # set modulation frequency
            self.frequency_set(freq_val_mhz * 1e6)

            # set up loop variables
            counter = 0
            timestamp_mu_list = [0] * self.num_counts
            self.core.break_realtime()

            # trigger sequence off same phase of RF
            self.rf_clock._set_sensitivity(1)
            time_trigger_rf_mu = self.rf_clock.timestamp_mu(now_mu() + self.time_rf_gating_mu)

            # start photon correlation sequence
            if time_trigger_rf_mu > 0:

                # set rtio hardware time to rising edge of RF
                at_mu(time_trigger_rf_mu + self.time_rf_holdoff_mu)
                self.rf_clock._set_sensitivity(0)

                # activate modulation
                at_mu(time_trigger_rf_mu + self.time_rf_holdoff_mu + self.time_rf_holdoff_mu)
                with parallel:
                    self.mod_toggle.on()
                    self.pmt_counter._set_sensitivity(1)
                    time_start_mu = now_mu() + self.time_mod_delay_mu

                #
                # # enable photon counting
                # at_mu(time_start_mu)
                # start counting photons
                while counter < self.num_counts:

                    # get photon timestamp
                    time_photon_mu = self.pmt_counter.timestamp_mu(now_mu() + self.time_pmt_gating_mu)

                    # move timestamped photon into buffer if valid
                    if time_photon_mu >= 0:
                        timestamp_mu_list[counter] = time_photon_mu
                        counter += 1

                # stop counting and upload
                self.core.break_realtime()
                with parallel:
                    self.pmt_counter._set_sensitivity(0)
                    self.mod_toggle.off()
                    self.update_dataset(time_start_mu, timestamp_mu_list)

            # if we don't get rf trigger for some reason, just reset
            else:
                self.rf_clock._set_sensitivity(0)

            # reset FIFOs
            self.core.reset()


    @rpc(flags={"async"})
    def update_dataset(self, time_start_mu, timestamp_mu_list):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.mutate_dataset(
            'results',
            self._dataset_counter,
            np.array(self.core.mu_to_seconds(np.array(timestamp_mu_list) - time_start_mu))
        )
        self._dataset_counter += 1


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        # set desired voltgae
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        sleep(0.2)

        # wait until voltage updates
        voltage_get_v = self.dc.voltage(channel)
        while np.abs(voltage_set_v - voltage_get_v) > 0.05:
            sleep(0.2)
            voltage_get_v = self.dc.voltage(channel)

        # print current voltage for verification
        print('\tvoltage set: {}'.format(voltage_get_v))

    @rpc
    def frequency_set(self, freq_hz):
        """
        Set the RF to the desired frequency.
        """
        freq_set_hz = self.fg.frequency(freq_hz)
        print('\tfrequency set: {}'.format(freq_set_hz))

    def fg_write(self, msg):
        """
        Write a GPIB message to the function generator.
        """
        self.fg.gpib_write(msg)

    def analyze(self):
        # turn off modulation
        self.fg.toggle(0)

        # process data
        # results_tmp = np.array(self.results).reshape((len(self.freq_mod_mhz_list), self.num_counts, 2))
        # ind_arr = np.argsort(self.freq_mod_mhz_list)
        # results_tmp = results_tmp[ind_arr]
        #
        # self.results_processed = results_tmp
