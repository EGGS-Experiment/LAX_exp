import numpy as np
from artiq.experiment import *

_DMA_HANDLE_INITIALIZE = "sideband_cooling_rabi_flop_initialize"
_DMA_HANDLE_SIDEBAND = "sideband_cooling_rabi_flop_pulse"
_DMA_HANDLE_READOUT = "sideband_cooling_rabi_flop_readout"


class SidebandCoolingRabiFlopping(EnvExperiment):
    """
    Sideband Cooling & Rabi Flopping
    Does sideband cooling then rabi flopping.
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
        self.setattr_argument("calibration",                                BooleanValue(default=False))
        self.setattr_argument("repetitions",                                NumberValue(default=25, ndecimals=0, step=1, min=1, max=10000))

        # additional cooling
        self.setattr_argument("repetitions_per_cooling",                    NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("additional_cooling_time_s",                  NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000))

        # sideband cooling
        self.setattr_argument("sideband_cycles",                            NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("cycles_per_spin_polarization",               NumberValue(default=150, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_min_sideband_cooling_us_list",          PYONValue([20]))
        self.setattr_argument("time_max_sideband_cooling_us_list",          PYONValue([200]))
        self.setattr_argument("freq_sideband_cooling_mhz_list",             PYONValue([104.118]))
        self.setattr_argument("ampl_sideband_cooling_pct",                  NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))

        # rabi-flopping readout
        self.setattr_argument("freq_qubit_mhz",                             NumberValue(default=104.463, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_rabi_us_list",                          Scannable(default=
                                                                                RangeScan(0, 200, 201, randomize=True),
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, ndecimals=5
                                                                            ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # ensure input has correct dimensions
        min_time_length = len(list(self.time_min_sideband_cooling_us_list))
        max_time_length = len(list(self.time_max_sideband_cooling_us_list))
        modes_length = len(list(self.freq_sideband_cooling_mhz_list))
        assert min_time_length == max_time_length == modes_length

        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_redist_mu =                                   self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_profileswitch_delay_mu =                      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_rfswitch_delay_mu =                           self.core.seconds_to_mu(2 * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # RF switches
        self.dds_repump_qubit_switch =                          self.get_device("ttl21")

        # convert frequency to ftw
        self.ftw_to_mhz =                                       1e3 / (2 ** 32 - 1)
        self.freq_redist_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                             self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                                   self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                             self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                   self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # sideband cooling
        self.time_sideband_cooling_list_mu =                    np.array([
                                                                    self.core.seconds_to_mu(time_us * us)
                                                                    for time_us in np.linspace(
                                                                        self.time_min_sideband_cooling_us_list,
                                                                        self.time_max_sideband_cooling_us_list,
                                                                        self.sideband_cycles
                                                                    )
                                                                ])

        # calculate number of spin polarizations
        num_spin_depolarizations = int(self.sideband_cycles / self.cycles_per_spin_polarization)
        if (num_spin_depolarizations < 1):
            num_spin_depolarizations = 1

        self.time_sideband_cooling_list_mu =                    np.array_split(self.time_sideband_cooling_list_mu, num_spin_depolarizations)

        # other sideband cooling parameters
        self.freq_sideband_cooling_ftw_list =                   [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_sideband_cooling_mhz_list]
        self.ampl_sideband_cooling_asf =                        self.dds_qubit.amplitude_to_asf(self.ampl_sideband_cooling_pct / 100)
        self.iter_sideband_cooling_modes_list =                 list(range(1, 1 + len(self.freq_sideband_cooling_ftw_list)))

        # rabi flopping timing
        max_time_us =                                           np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                            self.core.seconds_to_mu(np.array([
                                                                    [(max_time_us - time_us) * us, time_us * us]
                                                                    for time_us in self.time_rabi_us_list
                                                                ]))

        # calibration setup
        self.calibration_qubit_status =                         not self.calibration

        # set up datasets
        self.set_dataset("sideband_cooling_rabi_flopping", [])
        self.setattr_dataset("sideband_cooling_rabi_flopping")
        self.set_dataset("sideband_cooling_rabi_flopping_processed", np.zeros([len(self.time_rabiflop_mu_list), 3]))
        self.setattr_dataset("sideband_cooling_rabi_flopping_processed")


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        handle_initialize = self.core_dma.get_handle(_DMA_HANDLE_INITIALIZE)
        handle_sideband = self.core_dma.get_handle(_DMA_HANDLE_SIDEBAND)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep pi-pulse time
            for time_rabiflop_mu_list in self.time_rabiflop_mu_list:

                # get timings
                time_delay_mu = time_rabiflop_mu_list[0]
                time_rabi_mu = time_rabiflop_mu_list[1]

                # initialize by running doppler cooling and spin polarization
                self.core_dma.playback_handle(handle_initialize)

                # run sideband cooling cycles and repump afterwards
                self.core_dma.playback_handle(handle_sideband)

                # wait given time, and meanwhile set rabi-flopping waveform
                with parallel:
                    self.dds_qubit_board.set_profile(0)
                    delay_mu(time_delay_mu)

                # do rabi flopping
                self.dds_qubit.cfg_sw(True)
                delay_mu(time_rabi_mu)
                self.dds_qubit.cfg_sw(False)

                # read out
                self.core_dma.playback_handle(handle_readout)

                # record data
                self.update_dataset(time_rabi_mu, self.pmt_counter.fetch_count())
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
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)
        self.dds_repump_qubit_switch.on()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # doppler cooling sequence
        with self.core_dma.record(_DMA_HANDLE_INITIALIZE):

            # enable 854 rf switch
            # with parallel:
            self.dds_repump_qubit_switch.on()
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
                self.dds_qubit_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # doppler cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_doppler_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

        # sideband cooling sequence
        with self.core_dma.record(_DMA_HANDLE_SIDEBAND):

            # ensure state preparation is properly interspersed
            for time_list_mu in self.time_sideband_cooling_list_mu:

                # spin polarization/redistribute S-1/2 (397)
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

                # sweep pi-pulse times
                for time_modes_mu in time_list_mu:

                    # sweep over modes
                    for i in self.iter_sideband_cooling_modes_list:

                        # set mode
                        with parallel:
                            self.dds_qubit_board.set_profile(i)
                            delay_mu(self.time_profileswitch_delay_mu)

                        # qubit pi-pulse
                        self.dds_qubit.cfg_sw(self.calibration_qubit_status)
                        delay_mu(time_modes_mu[i - 1])
                        self.dds_qubit.cfg_sw(False)

                        # qubit repump
                        self.dds_board.cfg_switches(0b1100)
                        delay_mu(self.time_repump_qubit_mu)
                        self.dds_board.cfg_switches(0b0100)

            # repump qubit after sideband cooling
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):

            # set pump readout waveform
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

        # set AOM DDS waveforms
        # profile 0 = cooling; profile 1 = readout (red-detuned); profile 2 = readout (blue-detuned)
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

        # set rf switches
        self.dds_repump_qubit_switch.off()

        # set sideband cooling profiles
        # profile 0 = readout pi-pulse
        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=0)
        # profile 1 & greater = sideband cooling
        for i in self.iter_sideband_cooling_modes_list:
            self.dds_qubit.set_mu(self.freq_sideband_cooling_ftw_list[i - 1], asf=self.ampl_sideband_cooling_asf, profile=i)
            self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, time_mu, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('sideband_cooling_rabi_flopping', [self.core.mu_to_seconds(time_mu), pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # # turn dataset into numpy array for ease of use
        # self.sideband_cooling = np.array(self.sideband_cooling)
        #
        # # get sorted x-values (frequency)
        # freq_list_mhz = sorted(set(self.sideband_cooling[:, 0]))
        #
        # # collate results
        # collated_results = {
        #     freq: []
        #     for freq in freq_list_mhz
        # }
        # for freq_mhz, pmt_counts in self.sideband_cooling:
        #     collated_results[freq_mhz].append(pmt_counts)
        #
        # # process counts for mean and std and put into processed dataset
        # for i, (freq_mhz, count_list) in enumerate(collated_results.items()):
        #     binned_count_list = np.heaviside(np.array(count_list) - self.pmt_discrimination, 1)
        #     self.sideband_cooling_processed[i] = np.array([freq_mhz, np.mean(binned_count_list), np.std(binned_count_list)])
