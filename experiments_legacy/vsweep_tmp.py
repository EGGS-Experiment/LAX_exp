import labrad
import numpy as np
from time import sleep
from scipy import stats
from scipy import optimize

from os import environ
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE

from EGGS_labrad.config.dc_config import dc_config


class vsweeptmp(EnvExperiment):
    """
    v Sweep tmp
    """
    kernel_invariants = {
        'time_pmt_gating_mu',
        'dc_micromotion_channel',
        'dc_micromotion_voltage_v'
    }


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # num_counts
        self.setattr_argument("num_counts",                         NumberValue(default=10000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=22.5, ndecimals=1, step=0.5, min=0, max=31.5), group='mod')
        self.setattr_argument("mod_freq_khz",                       NumberValue(default=1701, ndecimals=4, step=1, min=0, max=400000), group='mod')


        # voltage
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel",             EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'), group='voltage')
        self.setattr_argument("dc_micromotion_voltages_v_list",     Scannable(
                                                                        default=CenterScan(60.0, 40.0, 1.0, randomize=True),
                                                                        global_min=0, global_max=400, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ), group='voltage')
        # self.setattr_argument("dc_micromotion_voltages_v_list",     Scannable(
        #                                                                 default=ExplicitScan([66]),
        #                                                                 global_min=0, global_max=400, global_step=1,
        #                                                                 unit="V", scale=1, ndecimals=4
        #                                                             ), group='voltage')

        # cooling
        self.setattr_argument("ampl_cooling_pct",                   NumberValue(default=35, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",                   NumberValue(default=100, ndecimals=6, step=1, min=1, max=500), group='cooling')


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl0")
        self.time_pmt_gating_mu =                                   self.core.seconds_to_mu(100 * us)
        self.pmt_flipper =                                          self.get_device("ttl23")

        # get voltage parameters
        self.dc_micromotion_channel_num =                           self.dc_micromotion_channeldict[self.dc_micromotion_channel]['num']
        self.dc_micromotion_channel_name =                          self.dc_micromotion_channel
        self.dc_micromotion_voltages_v_list =                       np.array(list(self.dc_micromotion_voltages_v_list))

        # cooling beam
        self.cooling_dds =                                          self.get_device("urukul1_ch1")
        self.cooling_dds_ampl_asf =                                 self.cooling_dds.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.cooling_dds_freq_ftw =                                 self.cooling_dds.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.cooling_dds_att_mu =                                   self.cooling_dds.cpld.att_to_mu(14 * dB)
        # cooling holdoff time
        self.time_cooling_holdoff_mu =                              self.core.seconds_to_mu(3 * ms)

        # modulation control and synchronization
        self.mod_dds =                                              self.get_device("urukul0_ch2")
        self.mod_dds_ampl_pct =                                     self.mod_dds.amplitude_to_asf(0.35)
        self.mod_dds_att_mu =                                       self.mod_dds.cpld.att_to_mu(self.mod_att_db * dB)
        self.mod_freq_mu =                                          self.mod_dds.frequency_to_ftw(self.mod_freq_khz * kHz)

        # RF synchronization
        self.rf_clock =                                             self.get_device('ttl7')
        self.time_rf_holdoff_mu =                                   self.core.seconds_to_mu(100000 * ns)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(150 * ns)


        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([1 * len(self.dc_micromotion_voltages_v_list),
                                                                              4]))
        self.setattr_dataset("results")

        # record parameters
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_attenuation_db',               self.mod_att_db)
        self.set_dataset('dc_channel_num',                          self.dc_micromotion_channel_num)
        self.set_dataset('dc_channel_name',                         self.dc_micromotion_channel_name)
        self.set_dataset('cooling_freq_mhz',                        self.freq_cooling_mhz)
        self.set_dataset('cooling_ampl_pct',                        self.ampl_cooling_pct)

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
        # sweep voltage
        for voltage_v in self.dc_micromotion_voltages_v_list:

            # reset timestamping loop counter
            counter = 0

            # set DC voltage
            self.voltage_set(self.dc_micromotion_channel_num, voltage_v)
            self.core.break_realtime()

            # add holdoff period for recooling the ion
            delay_mu(self.time_cooling_holdoff_mu)

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
                    self.update_dataset(self.mod_freq_mu, voltage_v, time_start_mu, timestamp_mu_list)

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
            self.pmt_flipper.off()

        # configure cooling dds
        self.cooling_dds.set_mu(self.cooling_dds_freq_ftw, asf=self.cooling_dds_ampl_asf)
        # self.cooling_dds.set_att_mu(self.cooling_dds_att_mu)
        self.cooling_dds.cpld.cfg_switches(0b1110)

        # configure rf modulation source
        self.mod_dds.cfg_sw(False)
        self.mod_dds.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.mod_dds.set_att_mu(self.mod_dds_att_mu)

        # set modulation frequency
        self.mod_dds.set_mu(self.mod_freq_mu, asf=self.mod_dds_ampl_pct)
        self.mod_dds.set_cfr1(phase_autoclear=1)

    @rpc
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)


    @rpc(flags={"async"})
    def update_dataset(self, freq_mu, voltage_v, time_start_mu, timestamp_mu_list):
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
            np.array([freq_mhz, voltage_v, np.abs(counts_demod), np.angle(counts_demod)])
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


    # ANALYSIS
    @rpc
    def _complexFitMinimize(self, dataset):
        """
        Extract the optimal voltage to minimize complex displacement
        from the RF origin.

        Arguments:
            dataset (list(list(float, complex)): the dataset comprised of the real
                independent variable, and the complex dependent variable.

        Returns:
            (float) : the value of the independent variable that minimizes the complex amplitude.
        """

        # create norm function for least squares optimization
        def func_norm(b_params, x, y):
            return b_params[0] + b_params[1] * x[1] - y

        # guess starting b_params
        b_guess = [0.1, 0.1]

        # do a complex least squares fit
        res = optimize.least_squares(func_norm, b_guess, args=(dataset[:, 0], dataset[:, 1]))
        res_intercept, res_slope = res.x

        # extract optimal voltage to minimize displacement
        voltage_optimal = - (np.re(res_intercept) * np.re(res_slope) + np.imag(res_intercept) * np.imag(res_slope)) / (
            np.power(np.abs(res_slope), 2))

        return voltage_optimal

    def analyze(self):
        """
        todo: document
        """
        # get results and format for complex linear fitting
        self.tmp0 = np.array(self.results, dtype='complex128')
        self.tmp0 = np.array([
            self.tmp0[:, 1],
            self.tmp0[:, 2] * np.exp(1.j * self.tmp0[:, 3])
        ], dtype='complex')
        self.tmp0 = self.tmp0.transpose()

        # extract minimum voltage
        min_voltage = self._complexFitMinimize(self.tmp0)

        # tmp remove
        print('\t\tmin voltage: {:f} V'.format(min_voltage))
