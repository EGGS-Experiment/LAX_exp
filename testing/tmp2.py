import numpy as np
from artiq.experiment import *
_DMA_HANDLE_LASERSCAN = "laserscan_sequence"


class LaserScanTrigger(EnvExperiment):
    """
    729nm Laser Scan Trigger
    Gets the number of counts as a function of frequency for a fixed time.
    """

    #kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",
        "time_profileswitch_delay_us",

        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",

        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_pump_rescue_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_pump_rescue_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct"
    ]


    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                            NumberValue(default=200, ndecimals=0, step=1, min=1, max=10000))

        # additional cooling
        self.setattr_argument("repetitions_per_cooling",                NumberValue(default=2000000, ndecimals=0, step=1, min=1, max=1000000000))
        self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000))

        # adc recording
        self.setattr_argument("adc_channel_gain_dict",                  PYONValue({0: 1000, 1: 10}))

        # timing
        self.setattr_argument("time_729_us",                            NumberValue(default=100, ndecimals=5, step=1, min=1, max=10000000))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",                    Scannable(
                                                                            default=CenterScan(103.77, 0.001, 0.001, randomize=True),
                                                                            global_min=60, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ))
        self.setattr_argument("att_729_db",                             NumberValue(default=22, ndecimals=1, step=0.5, min=8, max=32))


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # ADC
        self.adc =                                              self.get_device("sampler0")
        self.adc_channel_list =                                 list(self.adc_channel_gain_dict.keys())
        self.adc_gain_list_mu =                                 [int(np.log10(gain_mu)) for gain_mu in self.adc_channel_gain_dict.values()]
        self.adc_mu_to_v_list =                                 np.array([10 / (2**15 * gain_mu) for gain_mu in self.adc_channel_gain_dict.values()])

        # convert time values to machine units
        self.time_doppler_cooling_mu =                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_redist_mu =                                   self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_729_mu =                                      self.core.seconds_to_mu(self.time_729_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_profileswitch_delay_mu =                      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # process scan frequencies
        self.ftw_to_mhz =                                       1e3 / (2 ** 32 - 1)
        self.freq_qubit_scan_ftw =                              [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz]

        # convert dds values to machine units - frequency
        self.freq_redist_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                             self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert dds values to machine units - amplitude
        self.ampl_redist_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                             self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                   self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # attenuation
        self.att_729_mu =                                       self.dds_qubit.cpld.att_to_mu(self.att_729_db * dB)

        # set up datasets
        self._iter_dataset = 0
        self.set_dataset("laser_scan",                          np.zeros((self.repetitions * len(self.freq_qubit_scan_ftw), 2)))
        self.setattr_dataset("laser_scan")
        self.set_dataset("laser_scan_processed",                np.zeros([len(self.freq_qubit_scan_ftw), 3]))
        self.setattr_dataset("laser_scan_processed")
        self.set_dataset("adc_values",                          np.zeros((self.repetitions * len(self.freq_qubit_scan_ftw), 1 + len(self.adc_channel_list))))
        self.setattr_dataset("adc_values")

        # tmp remove
        self.setattr_device('ttl15')
        self.setattr_device('ttl22')
        print(list(self.freq_qubit_scan_mhz))


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # create ADC holding buffer
        sampler_buffer = [0] * 8
        self.core.break_realtime()

        # record dma
        self.DMArecord()
        self.core.break_realtime()

        # get dma handle
        handle = self.core_dma.get_handle(_DMA_HANDLE_LASERSCAN)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # set frequencies
            for freq_ftw in self.freq_qubit_scan_ftw:
                self.core.break_realtime()

                # set waveform for qubit DDS
                self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf)
                self.core.break_realtime()

                # run sequence
                self.core_dma.playback_handle(handle)

                # get ADC values
                self.adc.sample_mu(sampler_buffer)
                self.core.break_realtime()

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, self.pmt_counter.fetch_count(), now_mu(), sampler_buffer)
                    self.core.break_realtime()

            # add post repetition cooling
            if (trial_num > 0) and (trial_num % self.repetitions_per_cooling == 0):
                # set rescue waveform
                with parallel:
                    self.dds_board.set_profile(2)
                    delay_mu(self.time_profileswitch_delay_mu)

                # start rescuing
                self.dds_board.io_update.pulse_mu(8)
                self.dds_board.cfg_switches(0b0110)
                delay(self.additional_cooling_time_s)
                self.dds_board.cfg_switches(0b0100)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)

        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_board.io_update.pulse_mu(8)

        self.dds_qubit_board.set_profile(0)
        self.dds_qubit_board.io_update.pulse_mu(8)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE_LASERSCAN):

            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # repump pulse
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

            # cooling pulse
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_doppler_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # state preparation
            self.dds_board.cfg_switches(0b0101)
            delay_mu(self.time_redist_mu)
            self.dds_board.cfg_switches(0b0100)

            # 729
            self.ttl15.on()
            self.ttl22.off()
            self.dds_qubit.cfg_sw(True)
            delay_mu(self.time_729_mu)
            self.dds_qubit.cfg_sw(False)
            self.ttl15.off()
            self.ttl22.off()

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # readout pulse
            self.dds_board.cfg_switches(0b0110)
            self.pmt_gating_edge(self.time_readout_mu)
            self.dds_board.cfg_switches(0b0100)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # adjust attenuations correctly
        self.dds_qubit_board.get_att_mu()
        self.core.break_realtime()
        self.dds_qubit.set_att_mu(self.att_729_mu)
        self.core.break_realtime()

        # set AOM DDS waveforms
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=0)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=1)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=2)
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw, asf=self.ampl_pump_rescue_asf, profile=2)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=2)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=2)
        self.core.break_realtime()

        # set ADC channel gains
        for i in range(len(self.adc_channel_list)):
            self.adc.set_gain_mu(self.adc_channel_list[i], self.adc_gain_list_mu[i])
            self.core.break_realtime()

        # tmp remove
        self.dds_qubit.set_cfr1(ram_enable=0)
        self.dds_qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts, global_time_mu, adc_values_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # format adc data; convert values from mu to volts
        adc_values_volts = np.array(adc_values_mu)[self.adc_channel_list] * self.adc_mu_to_v_list
        adc_data = np.append(self.core.mu_to_seconds(global_time_mu), adc_values_volts)

        # save data to datasets
        self.mutate_dataset('laser_scan', self._iter_dataset, np.array([freq_ftw * self.ftw_to_mhz, pmt_counts]))
        self.mutate_dataset('adc_values', self._iter_dataset, adc_data)

        # update dataset iterator
        self._iter_dataset += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        for i in range(len(self.adc_channel_list)):
            print('\tch {:d}: {:.3f} +/- {:.3f} mV'.format(self.adc_channel_list[i], np.mean(self.adc_values[:, i + 1]) * 1000, np.std(self.adc_values[:, i + 1]) * 1000))
        # # tmp remove
        # self.pmt_discrimination = 17
        #
        # # turn dataset into numpy array for ease of use
        # self.laser_scan = np.array(self.laser_scan)
        #
        # # get sorted x-values (frequency)
        # freq_list_mhz = sorted(set(self.laser_scan[:, 0]))
        #
        # # collate results
        # collated_results = {
        #     freq: []
        #     for freq in freq_list_mhz
        # }
        # for freq_mhz, pmt_counts in self.laser_scan:
        #     collated_results[freq_mhz].append(pmt_counts)
        #
        # # process counts for mean and std and put into processed dataset
        # for i, (freq_mhz, count_list) in enumerate(collated_results.items()):
        #     binned_count_list = np.heaviside(np.array(count_list) - self.pmt_discrimination, 1)
        #     self.laser_scan_processed[i] = np.array([freq_mhz, np.mean(binned_count_list), np.std(binned_count_list)])