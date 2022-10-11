import numpy as np
from artiq.experiment import *

_DMA_HANDLE_DOPPLER = "sideband_cooling_doppler"
_DMA_HANDLE_SIDEBAND = "sideband_cooling_pulse"
_DMA_HANDLE_READOUT_RED = "sideband_cooling_readout_red"
_DMA_HANDLE_READOUT_BLUE = "sideband_cooling_readout_blue"


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
        "time_probe_us",
        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",
        "freq_probe_mhz",
        "freq_pump_cooling_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",
        "freq_qubit_mhz",
        "ampl_probe_pct",
        "ampl_pump_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",
        "att_probe_dB",
        "att_pump_dB",
        "att_pump_cooling_dB",
        "att_pump_readout_dB"
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
        self.setattr_argument("time_repetition_delay_us",       NumberValue(default=1000, ndecimals=5, step=1, min=1, max=10000))

        # AOM DDS parameters
        self.setattr_argument("freq_pump_readout_red_mhz",      NumberValue(default=92, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_pump_readout_blue_mhz",     NumberValue(default=92, ndecimals=3, step=1, min=10, max=200))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                  self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =              getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =      self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_sideband_cooling_mu =     self.core.seconds_to_mu(self.time_sideband_cooling_us * us)
        self.time_repump_qubit_mu =         self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_redist_mu =               self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_readout_mu =              self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_repetition_delay_mu =     self.core.seconds_to_mu(self.time_repetition_delay_us * us)

        # DDS devices
        self.dds_board =                    self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =              self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                    self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                     self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =           self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =             self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                    self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_probe_ftw =               self.dds_qubit.frequency_to_ftw(self.freq_probe_mhz * MHz)
        self.freq_pump_cooling_ftw =        self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_blue_ftw =   self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_blue_mhz * MHz)
        self.freq_pump_readout_red_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_red_mhz * MHz)
        self.freq_repump_cooling_ftw =      self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =        self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =               self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_probe_asf =               self.dds_qubit.amplitude_to_asf(self.ampl_probe_pct / 100)
        self.ampl_pump_asf =                self.dds_qubit.amplitude_to_asf(self.ampl_pump_pct / 100)
        self.ampl_repump_cooling_asf =      self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =        self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =               self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # sort out attenuation
        self.att_probe_mu =                 np.int32(0xFF) - np.int32(round(self.att_probe_dB * 8))
        self.att_pump_mu =                  np.int32(0xFF) - np.int32(round(self.att_pump_dB * 8))

        # set up datasets
        self.set_dataset("sideband_cooling", [])
        self.setattr_dataset("sideband_cooling")
        #self.set_dataset("parameters", [self.freq_qubit_ftw, self.att_cooling_mu, self.att_readout_mu])

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
        handle_doppler = self.core_dma.get_handle(_DMA_HANDLE_DOPPLER)
        handle_sideband = self.core_dma.get_handle(_DMA_HANDLE_SIDEBAND)
        handle_readout_red = self.core_dma.get_handle(_DMA_HANDLE_READOUT_RED)
        handle_readout_blue = self.core_dma.get_handle(_DMA_HANDLE_READOUT_BLUE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for i in range(self.repetitions):
            # run doppler cooling
            self.core_dma.playback_handle(handle_doppler)

            # run sideband cooling cycles
            for cycle_num in range(self.sideband_cycles):
                self.core_dma.playback_handle(handle_sideband)

            # do readout - red sideband
            self.core_dma.playback_handle(handle_readout_red)
            counts_red = self.pmt_counter.fetch_count()
            self.core.break_realtime()

            # do readout - blue sideband
            self.core_dma.playback_handle(handle_readout_blue)

            # record data & wait for ion to reheat
            with parallel:
                self.update_dataset(counts_red, self.pmt_counter.fetch_count())
                delay_mu(self.time_repetition_delay_mu)

        # reset DDSs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # doppler cooling sequence
        with self.core_dma.record(_DMA_HANDLE_DOPPLER):
            with sequential:
                # set cooling waveform
                self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_asf)

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

        # readout sequence - red
        with self.core_dma.record(_DMA_HANDLE_READOUT_RED):
            with sequential:
                # set readout waveform - red
                self.dds_pump.set_mu(self.freq_pump_readout_red_ftw, asf=self.ampl_pump_asf)

                # readout pulse
                self.dds_board.cfg_switches(0b0110)
                self.pmt_gating_edge(self.time_readout_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence - blue
        with self.core_dma.record(_DMA_HANDLE_READOUT_BLUE):
            with sequential:
                # set readout waveform - blue
                self.dds_pump.set_mu(self.freq_pump_readout_blue_ftw, asf=self.ampl_pump_asf)

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

        # set attenuations for now so we can change them later
        att_reg_old = np.int32(self.dds_board.get_att_mu())
        self.dds_board.set_all_att_mu(att_reg_old)
        self.core.break_realtime()

        # set AOM DDS waveforms
        self.dds_probe.set_att_mu(self.att_probe_mu)
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf)
        self.core.break_realtime()

        self.dds_pump.set_att_mu(self.att_pump_mu)
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_asf)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf)
        self.core.break_realtime()

        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, rsb_counts, bsb_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('sideband_cooling', [rsb_counts, bsb_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        self.sideband_cooling = np.array(self.sideband_cooling)
        print("red counts: {:.4f} +/- {:.4f}".format(np.mean(self.sideband_cooling[:, 0]), np.std(self.sideband_cooling[:, 0])))
        print("blue counts: {:.4f} +/- {:.4f}".format(np.mean(self.sideband_cooling[:, 1]), np.std(self.sideband_cooling[:, 1])))
