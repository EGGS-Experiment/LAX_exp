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
        'time_pmt_gating_mu',
        'dc_micromotion_channel_1',
        'dc_micromotion_channel_2',
        'dc_micromotion_voltages_v_list_1',
        'dc_micromotion_voltages_v_list_2',
        'mod_freq_mhz'
    }

    global_parameters = []


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # num_counts
        self.setattr_argument("num_counts",                         NumberValue(default=10000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=10, ndecimals=1, step=0.5, min=0, max=31.5))
        self.setattr_argument("mod_freq_mhz",                       NumberValue(default=1.415, ndecimals=5, step=0.001, min=0, max=1000000))


        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel_1",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_1",   Scannable(
                                                                        default=CenterScan(60.0, 40.0, 1.0, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        self.setattr_argument("dc_micromotion_channel_2",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='H Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_2",   Scannable(
                                                                        default=CenterScan(50.0, 40.0, 1.0, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl0")
        self.time_pmt_gating_mu =                                   self.core.seconds_to_mu(100 * us)

        # get voltage parameters
        self.dc_micromotion_voltages_v_list_1 =                     np.array(list(self.dc_micromotion_voltages_v_list_1))
        self.dc_micromotion_voltages_v_list_2 =                     np.array(list(self.dc_micromotion_voltages_v_list_2))
        self.dc_micromotion_channel_1_name =                        self.dc_micromotion_channel_1
        self.dc_micromotion_channel_1 =                             self.dc_micromotion_channeldict[self.dc_micromotion_channel_1]['num']
        self.dc_micromotion_channel_2_name =                        self.dc_micromotion_channel_2
        self.dc_micromotion_channel_2 =                             self.dc_micromotion_channeldict[self.dc_micromotion_channel_2]['num']

        # modulation control and synchronization
        self.mod_dds =                                              self.get_device("urukul0_ch2")
        self.mod_dds_ampl_pct =                                     self.mod_dds.amplitude_to_asf(0.35)
        self.mod_dds_att_mu =                                       self.mod_dds.cpld.att_to_mu(self.mod_att_db * dB)
        self.mod_freq_ftw =                                         self.mod_dds.frequency_to_ftw(self.mod_freq_mhz * MHz)


        # RF synchronization
        self.rf_clock =                                             self.get_device('ttl7')
        self.time_rf_holdoff_mu =                                   self.core.seconds_to_mu(100000 * ns)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(150 * ns)


        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([len(self.dc_micromotion_voltages_v_list_1) * len(self.dc_micromotion_voltages_v_list_2),
                                                                              4]))
        self.setattr_dataset("results")


        # record parameters
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_frequency_mhz',                self.mod_freq_mhz)
        self.set_dataset('modulation_attenuation_db',               self.mod_att_db)
        self.set_dataset('dc_channel_num_1',                        self.dc_micromotion_channel_1)
        self.set_dataset('dc_channel_num_2',                        self.dc_micromotion_channel_2)
        self.set_dataset('dc_channel_name_1',                       self.dc_micromotion_channel_1_name)
        self.set_dataset('dc_channel_name_2',                       self.dc_micromotion_channel_2_name)

        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                   self.cxn.dc_server


    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device
        self.core.reset()

        # prepare devices for experiment
        with parallel:
            self.prepareDevices()
            self.prepareDevicesLabrad()

        # set up loop variables
        counter = 0
        timestamp_mu_list = [0] * self.num_counts
        self.core.break_realtime()


        # MAIN LOOP
        # sweep voltage 1
        for voltage_1_v in self.dc_micromotion_voltages_v_list_1:

            # set voltage 1
            self.voltage_set(self.dc_micromotion_channel_1, voltage_1_v)
            self.core.break_realtime()

            # sweep voltage 2
            for voltage_2_v in self.dc_micromotion_voltages_v_list_2:

                # reset timestamping loop counter
                counter = 0

                # set voltage 2
                self.voltage_set(self.dc_micromotion_channel_2, voltage_2_v)
                self.core.break_realtime()

                # trigger sequence off same phase of RF
                self.rf_clock.gate_rising_mu(self.time_rf_gating_mu)
                time_trigger_rf_mu = self.rf_clock.timestamp_mu(now_mu())

                # start photon correlation sequence
                if time_trigger_rf_mu >= 0:

                    # activate modulation and enable photon counting
                    at_mu(time_trigger_rf_mu + self.time_rf_holdoff_mu)
                    self.mod_dds.cfg_sw(True)
                    with parallel:
                        self.pmt_counter._set_sensitivity(1)
                        time_start_mu = now_mu()
                        self.mod_dds.cpld.io_update.pulse_mu(8)

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
                        self.mod_dds.cfg_sw(False)
                        self.update_dataset(voltage_1_v, voltage_2_v, time_start_mu, timestamp_mu_list)

                # if we don't get rf trigger for some reason, just reset
                else:
                    self.rf_clock._set_sensitivity(0)

                # reset FIFOs
                self.core.reset()


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set ttl directions
        with parallel:
            self.pmt_counter.input()
            self.rf_clock.input()

        # configure rf modulation source
        self.mod_dds.cfg_sw(False)
        self.mod_dds.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.mod_dds.set_cfr1(phase_autoclear=1)
        self.core.break_realtime()

        # set rf modulation waveform
        self.mod_dds.set_mu(self.mod_freq_ftw, asf=self.mod_dds_ampl_pct)
        self.mod_dds.set_att_mu(self.mod_dds_att_mu)
        self.core.break_realtime()

    @rpc
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)


    @rpc(flags={"async"})
    def update_dataset(self, voltage_1_v, voltage_2_v, time_start_mu, timestamp_mu_list):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # remove starting time and digitally demodulate counts
        counts_mu = self.core.mu_to_seconds(np.array(timestamp_mu_list) - time_start_mu)
        counts_demod = np.sum(np.exp((2.j * np.pi * self.mod_freq_mhz * 1e6) * counts_mu)) / self.num_counts

        # update dataset
        self.mutate_dataset(
            'results',
            self._dataset_counter,
            np.array([voltage_1_v, voltage_2_v, np.abs(counts_demod), np.angle(counts_demod)])
        )
        self._dataset_counter += 1


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        voltage_set_v = self.dc.voltage_fast(channel, voltage_v)
        print('\tvoltage set: {}'.format(voltage_set_v))

    def analyze(self):
        pass
