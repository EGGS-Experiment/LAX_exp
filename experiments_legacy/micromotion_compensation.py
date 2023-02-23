import labrad
import numpy as np

from time import sleep
from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE
from EGGS_labrad.config.dc_config import dc_config


class MicromotionCompensation(EnvExperiment):
    """
    Micromotion Compensation

    Compensates micromotion by correlating photon counts with a modulation signal applied to the trapping RF.
    """
    kernel_invariants = {
        'time_timeout_pmt_mu',
        'dc_micromotion_channel_1',
        'dc_micromotion_channel_2',
        'dc_micromotion_voltages_v_list_1',
        'dc_micromotion_voltages_v_list_2',
        'ampl_mod_vpp',
        'freq_mod_mhz'
    }

    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge"
    ]


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # num_counts
        self.setattr_argument("num_counts",                         NumberValue(default=20000, ndecimals=0, step=1, min=1, max=10000000))

        # timing
        self.setattr_argument("time_timeout_pmt_s",                 NumberValue(default=50, ndecimals=5, step=1, min=1, max=1000000))

        # modulation
        self.setattr_argument("ampl_mod_vpp",                       NumberValue(default=0.5, ndecimals=3, step=0.01, min=0, max=1000000))
        self.setattr_argument("freq_mod_mhz",                       NumberValue(default=1.572, ndecimals=5, step=0.001, min=0, max=1000000))

        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel_1",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_1",   Scannable(
                                                                        default=CenterScan(70.0, 20.0, 0.1, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        self.setattr_argument("dc_micromotion_channel_2",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='H Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_2",   Scannable(
                                                                        default=CenterScan(45.0, 20.0, 0.1, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_timeout_pmt_mu =                                  self.core.seconds_to_mu(self.time_timeout_pmt_s * s)

        # get voltage parameters
        self.dc_micromotion_voltages_v_list_1 =                     np.array(list(self.dc_micromotion_voltages_v_list_1))
        self.dc_micromotion_voltages_v_list_2 =                     np.array(list(self.dc_micromotion_voltages_v_list_2))
        self.dc_micromotion_channel_1_name =                        self.dc_micromotion_channel_1
        self.dc_micromotion_channel_1 =                             self.dc_micromotion_channeldict[self.dc_micromotion_channel_1]['num']
        self.dc_micromotion_channel_2_name =                        self.dc_micromotion_channel_2
        self.dc_micromotion_channel_2 =                             self.dc_micromotion_channeldict[self.dc_micromotion_channel_2]['num']

        # RF modulation synchronization clock
        self.mod_clock =                                            self.get_device("urukul0_ch3")
        self.mod_toggle =                                           self.get_device("ttl8")
        self.mod_clock_freq_ftw =                                   self.mod_clock.frequency_to_ftw(10. * MHz)
        self.mod_clock_ampl_pct =                                   self.mod_clock.amplitude_to_asf(0.50)
        self.mod_clock_att_db =                                     0. * dB
        self.mod_clock_delay_mu =                                   self.core.seconds_to_mu(300 * ns)

        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([len(self.dc_micromotion_voltages_v_list_1) * len(self.dc_micromotion_voltages_v_list_2),
                                                                              self.num_counts]))
        self.setattr_dataset("results")


        # record parameters
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_frequency_hz',                 self.freq_mod_mhz)
        self.set_dataset('modulation_amplitude_vpp',                self.ampl_mod_vpp)
        self.set_dataset('channel_array',                           [self.dc_micromotion_channel_1_name, self.dc_micromotion_channel_2_name])
        self.set_dataset('voltage_array',                           np.stack(np.meshgrid(self.dc_micromotion_voltages_v_list_1,
                                                                                         self.dc_micromotion_voltages_v_list_2), -1).reshape(-1, 2))

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
        self.fg.frequency(self.freq_mod_mhz * 1e6)
        self.fg.burst(True)
        self.fg.burst_mode('GAT')
        self.fg.toggle(1)


    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # set ttl directions
        self.pmt_counter.input()
        self.mod_toggle.output()
        self.mod_toggle.off()
        self.core.break_realtime()

        # configure rf mod clock
        self.mod_clock.set_phase_mode(PHASE_MODE_ABSOLUTE)
        #self.core.break_realtime()
        self.mod_clock.set_att(self.mod_clock_att_db)
        self.mod_clock.set_mu(self.mod_clock_freq_ftw, asf=self.mod_clock_ampl_pct)
        self.mod_clock.cfg_sw(True)
        self.core.break_realtime()

        # enable external clocking
        self.fg_write(':ROSC:SOUR EXT')

        # MAIN LOOP
        # sweep voltage 1
        for voltage_1_v in self.dc_micromotion_voltages_v_list_1:

            # set voltage 1
            self.voltage_set(self.dc_micromotion_channel_1, voltage_1_v)
            self.core.break_realtime()

            # sweep voltage 2
            for voltage_2_v in self.dc_micromotion_voltages_v_list_2:

                # reset FIFOs
                self.core.reset()

                # set voltage 2
                self.voltage_set(self.dc_micromotion_channel_2, voltage_2_v)
                self.core.break_realtime()

                # set up loop variables
                counter = 0
                timestamp_mu_list = [0] * self.num_counts
                self.core.break_realtime()

                # synchronize timings with DDS clock
                #self.mod_clock.set_mu(self.mod_clock_freq_ftw, asf=self.mod_clock_ampl_pct)
                #delay_mu(self.mod_clock_delay_mu)

                # activate modulation and wait for change in output
                self.mod_toggle.on()
                delay_mu(self.mod_clock_delay_mu)
                time_start_mu = now_mu()

                # start counting photons
                time_stop_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
                while counter < self.num_counts:

                    # move timestamped photon into buffer
                    time_mu_tmp = self.pmt_counter.timestamp_mu(time_stop_mu)

                    # increase gating time
                    if time_mu_tmp < 0:
                        print('\t\tError: increased gating time.')
                        self.core.break_realtime()
                        time_stop_mu = self.pmt_counter.gate_rising_mu(self.time_timeout_pmt_mu)
                    else:
                        timestamp_mu_list[counter] = time_mu_tmp
                        counter += 1

                # stop loop
                self.pmt_counter._set_sensitivity(0)
                self.mod_toggle.off()
                self.core.reset()

                # store data
                self.update_dataset(time_start_mu, timestamp_mu_list)
                self.core.break_realtime()


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        # set desired voltgae
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        sleep(0.8)

        # wait until voltage updates
        voltage_get_v = self.dc.voltage(channel)
        while np.abs(voltage_set_v - voltage_get_v) > 0.05:
            sleep(0.8)
            print('voltage problem teehee')
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

    def analyze(self):
        # turn off modulation
        self.fg.toggle(0)

        # process data
        # results_tmp = np.array(self.results).reshape((len(self.dc_micromotion_voltages_v_list), self.num_counts, 2))
        # ind_arr = np.argsort(self.freq_mod_mhz_list)
        # results_tmp = results_tmp[ind_arr]
        #
        # self.results_processed = results_tmp
