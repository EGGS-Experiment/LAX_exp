import numpy as np
from random import shuffle
from artiq.experiment import *

_DMA_HANDLE_INITIALIZE = "eggs_heating_initialize"
_DMA_HANDLE_SIDEBAND = "eggs_heating_pulse"
_DMA_HANDLE_READOUT = "eggs_heating_readout"

# todo: make max time the same as readout time


class EGGSHeating(EnvExperiment):
    """
    EGGS Heating

    Does sideband cooling, then applies EGGS heating, then measures the sidebands.
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
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",

        "pmt_discrimination"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("calibration",                            BooleanValue(default=False))
        self.setattr_argument("repetitions",                            NumberValue(default=50, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("sideband_cycles",                        NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("cycles_per_spin_polarization",           NumberValue(default=150, ndecimals=0, step=1, min=1, max=10000))

        # sideband cooling
        self.setattr_argument("time_repump_sideband_cooling_us",        NumberValue(default=20, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_min_sideband_cooling_us_list",      PYONValue([20]))
        self.setattr_argument("time_max_sideband_cooling_us_list",      PYONValue([200]))
        self.setattr_argument("freq_sideband_cooling_mhz_list",         PYONValue([103.655]))
        self.setattr_argument("ampl_sideband_cooling_pct",              NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))

        # readout
        self.setattr_argument("freq_rsb_scan_mhz",                      Scannable(
                                                                            default=CenterScan(103.655, 0.04, 0.001),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ))

        self.setattr_argument("freq_bsb_scan_mhz",                      Scannable(
                                                                            default=CenterScan(105.271, 0.04, 0.001),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ))

        self.setattr_argument("time_readout_pipulse_us",                NumberValue(default=250, ndecimals=5, step=1, min=1, max=10000))

        # eggs heating
        self.setattr_argument("time_eggs_heating_ms",                   NumberValue(default=5, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
                                                                            default=CenterScan(85, 5, 0.2, randomize=True),
                                                                            global_min=30, global_max=400, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
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
        self.pmt_counter =                                              self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                                  self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_repump_qubit_mu =                                     self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_redist_mu =                                           self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                          self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_profileswitch_delay_mu =                              self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)

        # DDS devices
        self.dds_board =                                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                          self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.ftw_to_mhz =                                               1e3 / (2 ** 32 - 1)
        self.freq_redist_ftw =                                          self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # process scan frequencies
        self.freq_qubit_scan_ftw =                                      [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz)
                                                                        for freq_mhz in list(self.freq_rsb_scan_mhz) + list(self.freq_bsb_scan_mhz)]
        shuffle(self.freq_qubit_scan_ftw)

        # convert amplitude to asf
        self.ampl_redist_asf =                                          self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)

        # sideband cooling
        self.time_sideband_cooling_list_mu =                            np.array([
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

        self.time_sideband_cooling_list_mu =                            np.array_split(self.time_sideband_cooling_list_mu, num_spin_depolarizations)

        # other sideband cooling parameters
        self.time_repump_sideband_cooling_mu =                          self.core.seconds_to_mu(self.time_repump_sideband_cooling_us * us)
        self.freq_sideband_cooling_ftw_list =                           [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_sideband_cooling_mhz_list]
        self.ampl_sideband_cooling_asf =                                self.dds_qubit.amplitude_to_asf(self.ampl_sideband_cooling_pct / 100)
        self.iter_sideband_cooling_modes_list =                         list(range(1, 1 + len(self.freq_sideband_cooling_ftw_list)))

        # eggs heating
        self.awg_board =                                                self.get_device("phaser0")
        self.awg_eggs =                                                 self.awg_board.channel[0]
        self.time_eggs_heating_mu =                                     self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)
        # todo: get correct ftw conversion
        #self.awg_eggs_mhz_to_ftw =                                      lambda freq_mhz: np.int32(round(freq_mhz * ((1 << 30) / (6.25))))
        #self.freq_eggs_heating_ftw_list =                               [self.awg_eggs_mhz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_eggs_heating_mhz_list]

        # readout pi-pulse
        self.time_readout_pipulse_mu =                                  self.core.seconds_to_mu(self.time_readout_pipulse_us * us)
        self.ampl_qubit_asf =                                           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # calibration setup
        self.calibration_qubit_status =                                 not self.calibration

        # set up datasets
        self.set_dataset("eggs_heating", [])
        self.setattr_dataset("eggs_heating")
        self.set_dataset("eggs_heating_processed", np.zeros([len(self.freq_qubit_scan_ftw), 3]))
        self.setattr_dataset("eggs_heating_processed")


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
        for i in range(self.repetitions):

            # sweep eggs rf frequencies
            for freq_eggs_ftw in self.freq_eggs_heating_ftw_list:

                # todo: adjust dac nco freq
                # self.awg_eggs.oscillator[0].set_frequency(osc_freq_list[i] * MHz)

                # sweep final pi-pulse frequency
                for freq_ftw in self.freq_qubit_scan_ftw:

                    # set readout frequency in advance
                    self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf, profile=0)
                    self.core.break_realtime()

                    # initialize by running doppler cooling and spin polarization
                    self.core_dma.playback_handle(handle_initialize)

                    # run sideband cooling cycles and repump afterwards
                    self.core_dma.playback_handle(handle_sideband)

                    # pulse eggs heating
                    self.awg_eggs.oscillator[0].set_amplitude_phase(0.49, 0.)
                    self.awg_eggs.oscillator[1].set_amplitude_phase(0.49, 0.)
                    delay_mu(self.time_eggs_heating_mu)
                    self.awg_eggs.oscillator[0].set_amplitude_phase(0., 0.)
                    self.awg_eggs.oscillator[1].set_amplitude_phase(0., 0)

                    # read out
                    self.core_dma.playback_handle(handle_readout)

                    # record data
                    self.update_dataset(freq_ftw, self.pmt_counter.fetch_count(), freq_eggs_ftw)
                    self.core.break_realtime()

        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # doppler cooling sequence
        with self.core_dma.record(_DMA_HANDLE_INITIALIZE):

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
                        delay_mu(self.time_repump_sideband_cooling_mu)
                        self.dds_board.cfg_switches(0b0100)

            # repump qubit after sideband cooling
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

            # todo: do we need to state prep by spin polarization here?

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            # set pump readout waveform and qubit pi-pulse waveform
            with parallel:
                self.dds_board.set_profile(1)
                self.dds_qubit_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # do qubit pi-pulse
            self.dds_qubit.cfg_sw(True)
            delay_mu(self.time_readout_pipulse_mu)
            self.dds_qubit.cfg_sw(False)

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
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()

        # set sideband cooling profiles
        # profile 0 = readout pi-pulse, profile 1 & greater = sideband cooling
        for i in self.iter_sideband_cooling_modes_list:
            self.dds_qubit.set_mu(self.freq_sideband_cooling_ftw_list[i - 1], asf=self.ampl_sideband_cooling_asf, profile=i)
            self.core.break_realtime()


        # set up eggs awg (i.e. phaser)
        self.awg_board.init(debug=True)
        self.core.break_realtime()

        # nco stuff
        self.awg_eggs.set_nco_frequency(-212 * MHz)
        self.awg_eggs.set_nco_phase(0.)
        self.awg_board.dac_sync()
        self.core.break_realtime()

        # trf
        self.awg_eggs.set_att(0 * dB)
        self.awg_eggs.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()

        # duc
        self.awg_eggs.set_duc_frequency(0 * MHz)
        self.awg_eggs.set_duc_cfg()
        self.awg_board.duc_stb()
        self.core.break_realtime()

        # oscillators (i.e. sidebands)
        self.awg_eggs.oscillator[0].set_frequency(-self.freq_eggs_heating_secular_mhz * MHz)
        self.awg_eggs.oscillator[0].set_amplitude_phase(0.)

        self.awg_eggs.oscillator[1].set_frequency(self.freq_eggs_heating_secular_mhz * MHz)
        self.awg_eggs.oscillator[1].set_amplitude_phase(0.)


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts, time_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('eggs_heating', [freq_ftw * self.ftw_to_mhz, pmt_counts, self.core.mu_to_seconds(time_mu)])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.eggs_heating = np.array(self.eggs_heating)
        print('heating rate times')
        