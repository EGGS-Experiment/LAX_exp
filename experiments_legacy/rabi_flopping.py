import numpy as np
from artiq.experiment import *

_DMA_HANDLE_RESET = "rabi_flopping_reset"
_DMA_HANDLE_READOUT = "rabi_flopping_readout"


class RabiFlopping(EnvExperiment):
    """
    Rabi Flopping
    Measures ion fluorescence vs 729nm pulse time and frequency.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",

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
        self.setattr_argument("repetitions",                            NumberValue(default=40, ndecimals=0, step=1, min=1, max=10000))

        # additional cooling
        self.setattr_argument("repetitions_per_cooling",                NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000))

        # adc recording
        self.setattr_argument("adc_channel_gain_dict",                  PYONValue({0: 1000, 1: 10}))

        # qubit parameters
        self.setattr_argument("time_rabi_us_list",                      Scannable(
                                                                            default=RangeScan(0, 100, 401, randomize=True),
                                                                            global_min=1, global_max=100000, global_step=1,
                                                                            unit="us", scale=1, ndecimals=5
                                                                        ))

        # AOM values
        self.setattr_argument("freq_qubit_mhz",                         NumberValue(default=104.1055, ndecimals=5, step=1, min=1, max=10000))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                              self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # ADC
        self.adc =                                                      self.get_device("sampler0")
        self.adc_channel_list =                                         list(self.adc_channel_gain_dict.keys())
        self.adc_gain_list_mu =                                         [int(np.log10(gain_mu)) for gain_mu in self.adc_channel_gain_dict.values()]
        self.adc_mu_to_v_list =                                         np.array([10 / (2**15 * gain_mu) for gain_mu in self.adc_channel_gain_dict.values()])

        # convert time values to machine units
        self.time_profileswitch_delay_mu =                              self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_rfswitch_delay_mu =                                   self.core.seconds_to_mu(2 * us)
        self.time_repump_qubit_mu =                                     self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                                          self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_redist_mu =                                           self.core.seconds_to_mu(self.time_redist_us * us)

        # rabi flopping timing
        self.time_rabi_mu_list =                                        [self.core.seconds_to_mu(time_us * us) for time_us in self.time_rabi_us_list]
        max_time_us =                                                   np.max(list(self.time_rabi_us_list))
        self.time_delay_mu_list =                                       [self.core.seconds_to_mu((max_time_us - time_us) * us) for time_us in self.time_rabi_us_list]
        self.num_time_points_list =                                     list(range(len(self.time_rabi_mu_list)))

        # DDS devices
        self.dds_board =                                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                          self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # RF switches
        self.dds_repump_qubit_switch =                                  self.get_device("ttl21")

        # convert frequency to ftw
        self.freq_redist_ftw =                                          self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                                     self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                                           self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                          self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                                     self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # set up datasets
        self._iter_dataset =                                            0
        self.set_dataset("rabi_flopping",                               np.zeros((self.repetitions * len(self.time_rabi_mu_list), 2)))
        self.setattr_dataset("rabi_flopping")
        self.set_dataset("rabi_flopping_processed",                     np.zeros([len(self.time_rabi_mu_list), 2]))
        self.setattr_dataset("rabi_flopping_processed")
        self.set_dataset("adc_values",                                  np.zeros((self.repetitions * len(self.time_rabi_mu_list), 1 + len(self.adc_channel_list))))
        self.setattr_dataset("adc_values")


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

        # record dma and get handle
        self.DMArecord()
        handle_reset = self.core_dma.get_handle(_DMA_HANDLE_RESET)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep time
            for i in self.num_time_points_list:

                # get timings
                time_delay_mu = self.time_delay_mu_list[i]
                time_rabi_mu = self.time_rabi_mu_list[i]

                # run repump and cooling
                self.core_dma.playback_handle(handle_reset)

                # wait given time
                delay_mu(time_delay_mu)

                # rabi flopping w/qubit laser
                self.dds_qubit.cfg_sw(True)
                delay_mu(time_rabi_mu)
                self.dds_qubit.cfg_sw(False)

                # do readout
                self.core_dma.playback_handle(handle_readout)

                # get ADC values
                self.adc.sample_mu(sampler_buffer)
                self.core.break_realtime()

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_mu, self.pmt_counter.fetch_count(), now_mu(), sampler_buffer)
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

        # reset board profiles
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw, asf=self.ampl_pump_rescue_asf, profile=0)
        self.core.break_realtime()
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)
        self.dds_repump_qubit_switch.off()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # reset sequence
        with self.core_dma.record(_DMA_HANDLE_RESET):
            with sequential:
                # enable 854 rf switch
                # with parallel:
                self.dds_repump_qubit_switch.off()
                delay_mu(self.time_rfswitch_delay_mu)
                self.dds_board.cfg_switches(0b1100)

                # qubit repump (854) pulse
                delay_mu(self.time_repump_qubit_mu)
                # with parallel:
                self.dds_repump_qubit_switch.on()
                self.dds_board.cfg_switches(0b0100)
                delay_mu(self.time_rfswitch_delay_mu)

                # set cooling waveform
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # do cooling
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

                # do spin depolarization using probe
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
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

        # set AOM DDS waveforms; profile 0 is cooling, profile 1 is readout
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

        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf)
        self.core.break_realtime()

        # set rf switches
        self.dds_repump_qubit_switch.off()

        # set ADC channel gains
        for i in range(len(self.adc_channel_list)):
            self.adc.set_gain_mu(self.adc_channel_list[i], self.adc_gain_list_mu[i])
            self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, time_mu, pmt_counts, global_time_mu, adc_values_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        # format adc data; convert values from mu to volts
        adc_values_volts = np.array(adc_values_mu)[self.adc_channel_list] * self.adc_mu_to_v_list
        adc_data = np.append(self.core.mu_to_seconds(global_time_mu), adc_values_volts)

        # save data to datasets
        self.mutate_dataset('rabi_flopping', self._iter_dataset, np.array([self.core.mu_to_seconds(time_mu), pmt_counts]))
        self.mutate_dataset('adc_values', self._iter_dataset, adc_data)

        # update dataset iterator
        self._iter_dataset += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        for i in range(len(self.adc_channel_list)):
            print('\tch {:d}: {:.3f} +/- {:.3f} mV'.format(self.adc_channel_list[i],
                                                           np.mean(self.adc_values[:, i + 1]) * 1000,
                                                           np.std(self.adc_values[:, i + 1]) * 1000))
        # # tmp remove
        # self.pmt_discrimination = 17
        #
        # # turn dataset into numpy array for ease of use
        # self.rabi_flopping = np.array(self.rabi_flopping)
        #
        # # get sorted x-values (time, seconds)
        # time_list_s = sorted(set(self.rabi_flopping[:, 0]))
        #
        # # collate results
        # collated_results = {
        #     time: []
        #     for time in time_list_s
        # }
        # for time_s, pmt_counts in self.rabi_flopping:
        #     collated_results[time_s].append(pmt_counts)
        #
        # # process counts for mean and std and put into processed dataset
        # for i, (time_s, count_list) in enumerate(collated_results.items()):
        #     binned_count_list = np.heaviside(np.array(count_list) - self.pmt_discrimination, 1)
        #     self.rabi_flopping_processed[i] = np.array([time_s, np.mean(binned_count_list)])
