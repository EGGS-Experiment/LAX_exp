import labrad
import numpy as np
from time import sleep
from scipy import stats

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
        'time_pmt_gating_mu',
        'dc_micromotion_channel_1',
        'dc_micromotion_channel_2',
        'dc_micromotion_voltages_v_list_1',
        'dc_micromotion_voltages_v_list_2',
        'mod_freq_mhz'
    }


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # num_counts
        self.setattr_argument("num_counts",                         NumberValue(default=10000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=22.5, ndecimals=1, step=0.5, min=0, max=31.5), group='mod')
        self.setattr_argument("mod_freq_khz_list",                  Scannable(
                                                                        default=ExplicitScan([1.702, 1.545]),
                                                                        global_min=1, global_max=200000, global_step=1,
                                                                        unit="kHz", scale=1, ndecimals=4
                                                                    ), group='mod')
        self.setattr_argument("mod_freq_khz",                       NumberValue(default=1415, ndecimals=4, step=1, min=0, max=400000), group='mod')



        # voltage
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel_1",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_1",   Scannable(
                                                                        default=CenterScan(80.0, 60.0, 1.0, randomize=True),
                                                                        global_min=0, global_max=400, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        self.setattr_argument("dc_micromotion_channel_2",           EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='H Shim'))
        self.setattr_argument("dc_micromotion_voltages_v_list_2",   Scannable(
                                                                        default=CenterScan(80.0, 60.0, 1.0, randomize=True),
                                                                        global_min=0, global_max=400, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))

        # cooling
        self.setattr_argument("ampl_cooling_pct",                   NumberValue(default=35, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",                   NumberValue(default=100, ndecimals=6, step=1, min=1, max=500), group='cooling')


    def prepare(self):
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl0")
        self.time_pmt_gating_mu =                                   self.core.seconds_to_mu(100 * us)
        self.pmt_flipper =                                          self.get_device("ttl23")

        # get voltage parameters
        self.dc_micromotion_voltages_v_list_1 =                     np.array(list(self.dc_micromotion_voltages_v_list_1))
        self.dc_micromotion_voltages_v_list_2 =                     np.array(list(self.dc_micromotion_voltages_v_list_2))
        self.dc_micromotion_channel_1_name =                        self.dc_micromotion_channel_1
        self.dc_micromotion_channel_1_num =                         self.dc_micromotion_channeldict[self.dc_micromotion_channel_1]['num']
        self.dc_micromotion_channel_2_name =                        self.dc_micromotion_channel_2
        self.dc_micromotion_channel_2_num =                         self.dc_micromotion_channeldict[self.dc_micromotion_channel_2]['num']

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
        self.mod_freq_mu_list =                                     np.array([
                                                                        self.mod_dds.frequency_to_ftw(freq_mhz * kHz)
                                                                        for freq_mhz in self.mod_freq_khz_list
                                                                    ])

        # RF synchronization
        self.rf_clock =                                             self.get_device('ttl7')
        self.time_rf_holdoff_mu =                                   self.core.seconds_to_mu(100000 * ns)
        self.time_rf_gating_mu =                                    self.core.seconds_to_mu(150 * ns)


        # set up datasets
        self._dataset_counter                                       = 0
        self.set_dataset("results",                                 np.zeros([len(self.mod_freq_khz_list) * len(self.dc_micromotion_voltages_v_list),
                                                                              4]))
        self.setattr_dataset("results")

        # record parameters
        self.set_dataset('num_counts',                              self.num_counts)
        self.set_dataset('modulation_attenuation_db',               self.mod_att_db)
        self.set_dataset('dc_channel_1_num',                        self.dc_micromotion_channel_1_num)
        self.set_dataset('dc_channel_1_name',                       self.dc_micromotion_channel_1_name)
        self.set_dataset('dc_channel_2_num',                        self.dc_micromotion_channel_2_num)
        self.set_dataset('dc_channel_2_name',                       self.dc_micromotion_channel_2_name)
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


        # MAIN LOOP
        # sweep voltage
        for voltage_v in self.dc_micromotion_voltages_v_list:

            # set DC voltage
            self.voltage_set(self.dc_micromotion_channel_num, voltage_v)
            self.core.break_realtime()




    @kernel(flags={"fast-math"})
    def _sweep_voltage(self, voltage_freq_list):
        """
        todo: document
        Arguments:
            voltage_list_v  list(list(list(tuple(int, float)), float))): the list of voltages and frequencies to sweep.

        Returns:
                            list(list(float, float)): the correlated amplitude
                            and phase for each voltage.
        """
        # sweep over voltages
        for sweep_params in voltage_freq_list:

            # extract_params
            voltage_config_list = sweep_params[0]
            freq_mu = sweep_params[1]

            # prepare ions and devices
            with parallel:

                # set electrode voltages
                with sequential:

                    for voltage_params in voltage_config_list:
                        self.core.break_realtime()

                        # extract parameters
                        voltage_channel_num = voltage_params[0]
                        voltage_v = voltage_params[1]

                        # set voltage
                        self.voltage_set(self.dc_micromotion_channel_num, voltage_v)
                        self.core.break_realtime()

                # set modulation frequency
                with sequential:
                    self.mod_dds.set_mu(freq_mu, asf=self.mod_dds_ampl_pct)
                    self.mod_dds.set_cfr1(phase_autoclear=1)



            # add holdoff period for recooling the ion
            delay_mu(self.time_cooling_holdoff_mu)

            # get timestamped counts
            counts_timestamped_mu = self._timestamp_counts()

            # demodulate counts
            counts_demodulated_absarg = self._demodulate_counts(counts_timestamped_mu)


    @kernel(flags={"fast-math"})
    def _timestamp_counts(self):
        """
        ***todo
        """
        # set up loop variables
        counter = 0
        timestamp_mu_list = [0] * self.num_counts
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
            self.pmt_counter._set_sensitivity(0)
            self.mod_dds.cfg_sw(False)

        # if we don't get rf trigger for some reason, just reset
        else:
            self.rf_clock._set_sensitivity(0)

        # reset FIFOs
        self.core.reset()

        return timestamp_mu_list


    # PREPARE
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

    @rpc
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)


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
    @rpc(flags={"async"})
    def _demodulate_counts(self, freq_mu, time_start_mu, timestamp_mu_list):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # convert frequency to mhz
        freq_mhz = self.mod_dds.ftw_to_frequency(freq_mu) / MHz

        # remove starting time and digitally demodulate counts
        counts_mu = self.core.mu_to_seconds(np.array(timestamp_mu_list) - time_start_mu)
        counts_demod = np.sum(np.exp((2.j * np.pi * freq_mhz * 1e6) * counts_mu)) / len(timestamp_mu_list)

        # return demodulated counts
        return np.abs(counts_demod), np.angle(counts_demod)

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
        # split dataset into IV/DV, and real/imaginary
        dataset_x = dataset[:, 0]
        dataset_y = np.array([np.real(dataset[:, 1]), np.imag(dataset[:, 1])]).transpose()

        # fit the DV in the complex plane and get the unit vector of the line
        fit_complex = stats.linregress(dataset_y[:, 0], dataset_y[:, 1])
        m_c, b_c = fit_complex.slope, fit_complex.intercept

        # get the projection of the DV onto the fitted line
        vec_complex = np.linalg.norm(np.array([1, m_c]))
        dataset_proj = np.dot(dataset_y, vec_complex)

        # fit the x dataset to the parameterized line
        fit_parameterized = stats.linregress(dataset_x, dataset_proj)
        m_p, b_p = fit_complex.slope, fit_complex.intercept

        # extract x_min
        x_min = - (m_c * b_c)/(1 + np.pow(m_c, 2))

        # convert x_min to V_min
        V_min = (x_min - b_p) / m_p

        return V_min

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


    def analyze(self):
        pass
