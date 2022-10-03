import numpy as np
from artiq.experiment import *

# todo: turn relevant lasers off before switching and delay, then turn on
_DMA_HANDLE_COOLING = "rabi_flopping_cooling"
_DMA_HANDLE_READOUT = "rabi_flopping_readout"
_DMA_HANDLE_REPUMP = "rabi_flopping_repump"


class RabiFloppingRDX(EnvExperiment):
    """
    Rabi Flopping RDX
    Measures ion fluorescence vs 729nm pulse time and frequency, but better than before
    """

    # kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=2, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_profileswitch_delay_us",    NumberValue(default=500, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_readout_us",                NumberValue(default=100, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_repump_qubit_us",           NumberValue(default=100, ndecimals=5, step=1, min=1, max=10000))

        # AOM DDS channels
        self.setattr_argument("dds_board_num",                  NumberValue(default=1, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_probe_channel",              NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_pump_channel",               NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_repump_cooling_channel",     NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_repump_qubit_channel",       NumberValue(default=3, ndecimals=0, step=1, min=0, max=3))

        # AOM DDS channels - qubit
        self.setattr_argument("dds_board_qubit_num",            NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_qubit_channel",              NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

        # AOM DDS parameters
        self.setattr_argument("freq_probe_mhz",                 NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_pump_cooling_mhz",          NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_pump_readout_mhz",          NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_repump_cooling_mhz",        NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_repump_qubit_mhz",          NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_qubit_mhz",                 NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))

        self.setattr_argument("ampl_probe_pct",                 NumberValue(default=50, ndecimals=3, step=1, min=1, max=100))
        self.setattr_argument("ampl_pump_pct",                  NumberValue(default=50, ndecimals=3, step=1, min=1, max=100))
        self.setattr_argument("ampl_repump_cooling_pct",        NumberValue(default=50, ndecimals=3, step=1, min=1, max=100))
        self.setattr_argument("ampl_repump_qubit_pct",          NumberValue(default=50, ndecimals=3, step=1, min=1, max=100))
        self.setattr_argument("ampl_qubit_pct",                 NumberValue(default=50, ndecimals=3, step=1, min=1, max=100))

        self.setattr_argument("att_pump_cooling_dB",            NumberValue(default=21.5, ndecimals=1, step=0.5, min=8, max=31.5))
        self.setattr_argument("att_pump_readout_dB",            NumberValue(default=21.5, ndecimals=1, step=0.5, min=8, max=31.5))

        # qubit time scan
        self.setattr_argument("time_sweep_us",                  Scannable(default=RangeScan(1, 1000, 21),
                                                                    global_min=1, global_max=10000, global_step=1,
                                                                    unit="us", scale=1, ndecimals=0))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(default=RangeScan(109, 111, 401),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=3))

        # PMT
        self.setattr_argument("pmt_input_channel",              NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",                EnumerationValue(["rising", "falling", "both"], default="rising"))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_profileswitch_delay_mu =  self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_readout_mu =              self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_repump_qubit_mu =         self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_sweep_mu =                [self.core.seconds_to_mu(time_us * us) for time_us in self.time_sweep_us]

        # DDS devices
        self.dds_board =            self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =      self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =            self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =             self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =   self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =     self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =            self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert dds values to machine units - frequency
        self.ftw_to_mhz = 1e3 / (2 ** 32 - 1)
        #self.freq_qubit_scan_mu = list(self.freq_qubit_scan_mhz)
        # todo: check that list comprehension works
        self.freq_qubit_scan_ftw = [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz]

        self.freq_probe_ftw =           self.dds_qubit.frequency_to_ftw(self.freq_probe_mhz * MHz)
        self.freq_pump_cooling_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =           self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # convert dds values to machine units - amplitude
        self.ampl_probe_asf =           self.dds_qubit.amplitude_to_asf(self.ampl_probe_pct / 100)
        self.ampl_pump_asf =            self.dds_qubit.amplitude_to_asf(self.ampl_pump_pct / 100)
        self.ampl_repump_cooling_asf =  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # sort out attenuation
        self.att_cooling_reg = np.int32(0xFF)
        self.att_readout_reg = np.int32(0xFF)
        self.att_cooling_mu = np.int32(0xFF) - np.int32(round(self.att_pump_cooling_dB * 8))
        self.att_readout_mu = np.int32(0xFF) - np.int32(round(self.att_pump_readout_dB * 8))

        # set up datasets
        self.set_dataset("rabi_flopping", [], broadcast=True)

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
        handle = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # set frequencies
            for freq_ftw in self.freq_qubit_scan_ftw:
                self.core.break_realtime()

                # set waveform for qubit DDS
                self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf)
                self.core.break_realtime()

                # sweep time
                for rabi_time_mu in self.time_sweep_mu:
                    # turn on qubit repump (854) to pump back down
                    with parallel:
                        self.dds_board.cfg_switches(0b1100)
                        delay_mu(self.time_repump_qubit_mu)
                    self.dds_board.cfg_switches(0b0100)

                    # get pmt counts w/397 on for calibration
                    self.core_dma.playback_handle(handle)
                    pmt_calib = self.pmt_counter.fetch_count()
                    self.core.break_realtime()

                    # rabi flopping w/qubit laser
                    with parallel:
                        self.dds_qubit.cfg_sw(1)
                        delay_mu(rabi_time_mu)
                    self.dds_qubit.cfg_sw(0)

                    # get pmt counts (actual)
                    self.core_dma.playback_handle(handle)

                    # update dataset
                    with parallel:
                        self.update_dataset(rabi_time_mu, freq_ftw, self.pmt_counter.fetch_count(), pmt_calib)
                        self.core.break_realtime()

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
                # set cooling attenuation
                self.dds_pump.set_att_mu(self.att_cooling_mu)

                # change profile to allow readout light
                with parallel:
                    self.dds_board.set_profile(1)
                    self.dds_qubit_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # read PMT counts
                self.pmt_gating_edge(self.time_readout_mu)

                # set readout attenuation
                self.dds_pump.set_att_mu(self.att_readout_mu)

                # change profile to stop readout light
                with parallel:
                    self.dds_board.set_profile(0)
                    self.dds_qubit_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # initialize dds boards
        # self.dds_board.init()
        self.core.break_realtime()
        # self.dds_qubit_board.init()
        self.core.break_realtime()

        # set attenuations for now so we can change them later
        att_reg_old = np.int32(self.dds_board.get_att_mu())
        self.dds_board.set_all_att_mu(att_reg_old)
        self.core.break_realtime()

        # set AOM DDS waveforms
        # profile 0 is block, profile 1 is pass
        # self.dds_probe.init()
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf, profile=0)
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf, profile=1)
        self.dds_probe.cfg_sw(0)
        self.core.break_realtime()

        # self.dds_pump.init()
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_asf, profile=1)
        self.dds_pump.cfg_sw(1)
        self.core.break_realtime()

        # self.dds_repump_cooling.init()
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.dds_repump_cooling.cfg_sw(1)
        self.core.break_realtime()

        # self.dds_repump_qubit.init()
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.dds_repump_qubit.cfg_sw(1)
        self.core.break_realtime()

        # self.dds_qubit.init()
        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.dds_qubit.cfg_sw(0)
        self.core.break_realtime()

        # set correct profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)
        self.core.break_realtime()

    @rpc(flags={"async"})
    def update_dataset(self, time_mu, freq_ftw, pmt_counts, pmt_calib):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('rabi_flopping', [self.core.mu_to_seconds(time_mu), freq_ftw * self.ftw_to_mhz, pmt_counts, pmt_calib])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # todo: see if we can process data here instead, e.g. converting units
