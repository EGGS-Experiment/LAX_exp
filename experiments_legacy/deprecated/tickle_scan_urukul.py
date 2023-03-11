import labrad
import numpy as np
from os import environ
from time import sleep
from artiq.experiment import *

_DMA_HANDLE_TICKLE_URUKUL_SEQUENCE = "tickle_scan_urukul_sequence"


class TickleScanUrukul(EnvExperiment):
    """
    Tickle Scan Urukul
    Applies RF from an Urukul channel in search of secular frequencies and measures fluorescence.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "dds_board_tickle_num",
        "dds_tickle_channel",

        "freq_pump_cooling_mhz",
        "freq_repump_cooling_mhz",

        "ampl_pump_cooling_pct",
        "ampl_repump_cooling_pct"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                        NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # tickle values
        self.setattr_argument("time_tickle_us",                     NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("ampl_tickle_pct",                    NumberValue(default=50, ndecimals=3, step=1, min=10, max=10000))
        self.setattr_argument("freq_tickle_mhz",                    Scannable(
                                                                        default=RangeScan(0.9, 1.2, 31),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=4
                                                                    ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # DDS devices
        self.dds_board =                                            self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_pump =                                             self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                   self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))

        # convert frequency to ftw
        self.freq_pump_cooling_ftw =                                self.dds_pump.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_repump_cooling_ftw =                              self.dds_pump.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)

        # convert amplitude to asf
        self.ampl_pump_cooling_asf =                                self.dds_pump.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_repump_cooling_asf =                              self.dds_pump.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

        # tickle
        # self.dds_board_tickle =                                     self.get_device("urukul{:d}_cpld".format(self.dds_board_tickle_num))
        self.dds_board_tickle =                                     self.get_device("urukul{:d}_cpld".format(0))
        self.dds_tickle =                                           self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_tickle_num, self.dds_tickle_channel))
        self.time_tickle_mu =                                       self.core.seconds_to_mu(self.time_tickle_us * us)
        self.freq_tickle_ftw =                                      [self.dds_pump.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in list(self.freq_tickle_mhz)]
        self.ampl_tickle_asf =                                      self.dds_pump.amplitude_to_asf(self.ampl_tickle_pct / 100)

        # set up datasets
        self.set_dataset("tickle_scan_urukul", [])
        self.setattr_dataset("tickle_scan_urukul")


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
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_TICKLE_URUKUL_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # scan frequency
            for freq_ftw in self.freq_tickle_ftw:

                # set frequency
                self.dds_tickle.set_mu(freq_ftw, asf=self.ampl_tickle_asf)

                # run sequence
                self.core_dma.playback_handle(handle_sequence)

                # add data to dataset
                with parallel:
                    self.update_dataset(freq_ftw, self.pmt_counter.fetch_count())
                    self.core.break_realtime()

                # add artificial delay
                delay_mu(10 * ms)


        # reset devices
        self.dds_board.set_profile(0)

        self.dds_board.cfg_switches(0b1110)
        self.dds_board_tickle.cfg_switches(0b0000)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # tickle sequence
        with self.core_dma.record(_DMA_HANDLE_TICKLE_URUKUL_SEQUENCE):

            # start tickle
            self.dds_tickle.cfg_sw(True)

            # read pmt
            self.pmt_gating_edge(self.time_tickle_mu)

            # stop tickle
            self.dds_tickle.cfg_sw(False)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up AOMs
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('tickle_scan_urukul', [self.dds_tickle.ftw_to_frequency(freq_ftw), pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        self.tickle_scan_urukul = np.array(self.tickle_scan_urukul)