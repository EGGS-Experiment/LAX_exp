import numpy as np
from artiq.experiment import *

_DMA_HANDLE_INITIALIZE = "sideband_cooling_initialize"
_DMA_HANDLE_SIDEBAND = "sideband_cooling_pulse"
_DMA_HANDLE_READOUT = "sideband_cooling_readout"


class SidebandCooling(EnvExperiment):
    """
    Sideband Cooling
    Measures temperature after a given number of RSB pulses.
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
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("sideband_cycles",                NumberValue(default=4, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_sideband_cooling_us",       NumberValue(default=200, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_pipulse_us",                NumberValue(default=100, ndecimals=5, step=1, min=1, max=10000))

        # AOM
        self.setattr_argument("freq_sideband_cooling_mhz",      NumberValue(default=110, ndecimals=5, step=1, min=1, max=10000))

        self.setattr_argument("ampl_sideband_cooling_pct",      NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))
        self.setattr_argument("ampl_pipulse_pct",               NumberValue(default=50, ndecimals=5, step=1, min=1, max=100))


        # frequency sweep
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(default=
                                                                          RangeScan(104.24, 104.96, 801),
                                                                          global_min=30, global_max=200, global_step=1,
                                                                          unit="MHz", scale=1, ndecimals=5))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_redist_mu =                                   self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_redist_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)

        # novel parameters
        self.time_sideband_cooling_mu =                         self.core.seconds_to_mu(self.time_sideband_cooling_us * us)
        self.time_pipulse_mu =                                  self.core.seconds_to_mu(self.time_pipulse_us * us)

        self.freq_sideband_cooling_ftw =                        self.dds_qubit.frequency_to_ftw(self.freq_sideband_cooling_mhz * MHz)

        self.ampl_sideband_cooling_asf =                        self.dds_qubit.amplitude_to_asf(self.ampl_sideband_cooling_pct / 100)
        self.ampl_pipulse_asf =                                 self.dds_qubit.amplitude_to_asf(self.ampl_pipulse_pct / 100)

        # set up datasets
        self.set_dataset("sideband_cooling", [])
        self.setattr_dataset("sideband_cooling")


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

            # frequency sweep
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set readout frequency in advance
                self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_pipulse_asf, profile=1)
                self.core.break_realtime()

                # initialize by running doppler cooling and spin polarization
                self.core_dma.playback_handle(handle_initialize)

                # run sideband cooling cycles
                for cycle_num in range(self.sideband_cycles):
                    self.core_dma.playback_handle(handle_sideband)

                # repump qubit
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

                # read out
                self.core_dma.playback_handle(handle_readout)

                # record data
                self.update_dataset(freq_ftw, self.pmt_counter.fetch_count())
                self.core.break_realtime()

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # doppler cooling sequence
        with self.core_dma.record(_DMA_HANDLE_INITIALIZE):
            with sequential:
                # set qubit to sideband waveform
                self.dds_qubit_board.set_profile(0)

                # set cooling waveform
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # doppler cooling
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_doppler_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

                # spin depolarization/redistribute S-1/2 (397)
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

        # sideband cooling sequence
        with self.core_dma.record(_DMA_HANDLE_SIDEBAND):
            with sequential:
                # sideband cooling (729)
                self.dds_qubit.cfg_sw(1)
                delay_mu(self.time_sideband_cooling_mu)
                self.dds_qubit.cfg_sw(0)

                # qubit repump (854) pulse
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

                # spin depolarization/redistribute S-1/2 (397)
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
                # set qubit pi-pulse waveform
                self.dds_qubit_board.set_profile(1)

                # set pump readout waveform
                with parallel:
                    self.dds_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # do qubit pi-pulse
                self.dds_qubit.cfg_sw(1)
                delay_mu(self.time_pipulse_mu)
                self.dds_qubit.cfg_sw(0)

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
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()

        self.dds_qubit.set_mu(self.freq_sideband_cooling_mhz, asf=self.ampl_sideband_cooling_mhz, profile=0)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('sideband_cooling', [freq_ftw * self.ftw_to_frequency, pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        self.sideband_cooling = np.array(self.sideband_cooling)
