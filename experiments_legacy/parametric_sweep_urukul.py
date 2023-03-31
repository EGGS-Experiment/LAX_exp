import labrad
import numpy as np
from time import sleep

from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from EGGS_labrad.config.dc_config import dc_config


class ParametricSweepUrukul(EnvExperiment):
    """
    Parametric Sweep Urukul
    """
    kernel_invariants = {
        'time_pmt_gating_mu',
        'dc_micromotion_channel',
        'ampl_mod_vpp',
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
        self.setattr_argument("num_counts",                         NumberValue(default=10000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=30, ndecimals=1, step=0.5, min=5, max=31.5))
        self.setattr_argument("mod_freq_mhz_list",                  Scannable(
                                                                        default=CenterScan(1.209, 0.02, 0.0004, randomize=True),
                                                                        global_min=0, global_max=1000, global_step=0.001,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ))


        # voltage values
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel",             EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltage_v",           NumberValue(default=37.0, ndecimals=3, step=1, min=0, max=1000000))


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

        # modulation control and synchronization
        self.mod_dds =                                              self.get_device("urukul0_ch2")
        self.mod_dds_ampl_pct =                                     self.mod_dds.amplitude_to_asf(0.35)
        self.mod_dds_att_db =                                       10. * dB
        self.mod_freq_mu_list =                                     np.array([
                                                                        self.mod_dds.frequency_to_ftw(freq_mhz * MHz)
                                                                        for freq_mhz in self.mod_freq_mhz_list
                                                                    ])
        self.time_mod_delay_mu =                                    self.core.seconds_to_mu(1000 * ns)

        # RF synchronization
        self.rf_clock =                                             self.get_device('ttl7')
        self.time_rf_holdoff_mu =                                   self.core.seconds_to_mu(10000 * ns)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(100 * ns)

        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([len(self.mod_freq_mhz_list), 3]))
        self.setattr_dataset("results")

        # record parameters
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_attenuation_db',               self.mod_att_db)
        self.set_dataset('dc_channel_n',                            self.dc_micromotion_channel)
        self.set_dataset('dc_channel_voltage',                      self.dc_micromotion_voltage_v)

        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                   self.cxn.dc_server


    @kernel(flags={"fast-math"})
    def run(self):
        # reset core device
        self.core.reset()

        # prepare devices for experiment
        self.prepareDevices()

        # set up loop variables
        counter = 0
        timestamp_mu_list = [0] * self.num_counts
        self.core.break_realtime()


        # MAIN LOOP
        # sweep modulation frequency
        for freq_mu in self.mod_freq_mu_list:

            # reset timestamping loop counter
            counter = 0

            # set modulation frequency
            self.mod_dds.set_mu(freq_mu, asf=self.mod_dds_ampl_pct)
            self.mod_dds.set_cfr1(phase_autoclear=1)

            # trigger sequence off same phase of RF
            self.rf_clock._set_sensitivity(1)
            time_trigger_rf_mu = self.rf_clock.timestamp_mu(now_mu() + self.time_rf_gating_mu)

            # start photon correlation sequence
            if time_trigger_rf_mu >= 0:

                # set rtio hardware time to rising edge of RF
                at_mu(time_trigger_rf_mu + self.time_rf_holdoff_mu)
                self.rf_clock._set_sensitivity(0)

                # activate modulation and enable photon counting
                at_mu(time_trigger_rf_mu + 4 * self.time_rf_holdoff_mu)
                with parallel:
                    self.pmt_counter._set_sensitivity(1)
                    with sequential:
                        self.mod_dds.cfg_sw(True)
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
                # self.core.break_realtime()
                with parallel:
                    self.pmt_counter._set_sensitivity(0)
                    self.mod_dds.cfg_sw(False)
                    self.update_dataset(freq_mu, time_start_mu, timestamp_mu_list)

            # if we don't get rf trigger for some reason, just reset
            else:
                self.rf_clock._set_sensitivity(0)
                self.core.break_realtime()

            # reset FIFOs
            self.core.reset()


        # finish



    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set ttl directions
        with parallel:
            self.pmt_counter.input()
            self.rf_clock.input()
        self.core.break_realtime()

        # configure rf mod clock
        self.mod_dds.cfg_sw(False)
        self.mod_dds.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.mod_dds.set_att(self.mod_dds_att_db)
        self.core.break_realtime()

        # prepare LabRAD devices
        self._prepareDevicesLabrad()
        self.core.break_realtime()

    @rpc
    def _prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)

        # set voltage
        self.voltage_set(self.dc_micromotion_channel, self.dc_micromotion_voltage_v)


    @rpc(flags={"async"})
    def update_dataset(self, freq_mu, time_start_mu, timestamp_mu_list):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # convert frequency to mhz
        freq_mhz = self.mod_dds.ftw_to_frequency(freq_mu) / MHz

        # remove starting time and digitally demodulate counts
        counts_mu = self.core.mu_to_seconds(np.array(timestamp_mu_list) - time_start_mu)
        counts_demod = np.sum(np.exp((2.j * np.pi * freq_mhz * 1e6) * counts_mu)) / self.num_counts

        # update dataset
        self.mutate_dataset(
            'results',
            self._dataset_counter,
            np.array([freq_mhz, np.abs(counts_demod), np.angle(counts_demod)])
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
