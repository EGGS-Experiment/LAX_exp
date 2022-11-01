import labrad
import numpy as np
from os import environ
from time import sleep
from artiq.experiment import *

_DMA_HANDLE_TICKLE_SEQUENCE = "tickle_scan_sequence"


class TickleScan(EnvExperiment):
    """
    Tickle Scan
    Scans the RF from a function generator and measures fluorescence.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",
        "ttl_channel_function_generator",

        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "freq_pump_cooling_mhz",
        "freq_repump_cooling_mhz",

        "ampl_pump_cooling_pct",
        "ampl_repump_cooling_pct",

        "pmt_discrimination"
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
        self.setattr_argument("time_freq_delay_s",                  NumberValue(default=0.5, ndecimals=3, step=0.1, min=0, max=10))
        self.setattr_argument("ampl_tickle_mvpp",                   NumberValue(default=50, ndecimals=3, step=1, min=10, max=10000))
        self.setattr_argument("freq_tickle_mhz",                    Scannable(
                                                                        default=RangeScan(0.9, 1.2, 31),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=4
                                                                    ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.function_generator_server


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                          self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                      getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # function generator
        self.ttl_function_generator =                               self.get_device("ttl{}".format(self.ttl_channel_function_generator))

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
        self.time_fuction_generator_delay_mu =                      self.core.seconds_to_mu(350 * ns)
        self.time_tickle_us =                                       self.core.seconds_to_mu(self.time_tickle_us * us)
        self.freq_tickle_mhz =                                      list(self.freq_tickle_mhz)

        # prepare labrad devices
        self.fg.select_device()
        self.fg.toggle(1)
        self.fg.amplitude(self.ampl_tickle_mvpp / 1e3)

        # configure burst mode
        self.fg.gpib_write('BURS:STAT ON')
        self.fg.gpib_write('BURS:MODE GAT')
        self.fg.gpib_write('BURS:PHAS 0')

        # configure function generator trigger
        self.fg.gpib_write('TRIG:SOUR EXT')
        self.fg.gpib_write('TRIG:SLOP POS')

        # set up datasets
        self.set_dataset("tickle_scan", [])
        self.setattr_dataset("tickle_scan")


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
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_TICKLE_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        # for freq_mhz in self.freq_tickle_mhz:
        #
        #     # set frequency
        #     self.frequency_set(freq_mhz)
        #
        #     # repeat experiment
        #     for trial_num in range(self.repetitions):
        #
        #         # run sequence
        #         self.core_dma.playback_handle(handle_sequence)
        #
        #         # add data to dataset
        #         with parallel:
        #             self.update_dataset(freq_mhz, self.pmt_counter.fetch_count())
        #             self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # repeat experiment
            for freq_mhz in self.freq_tickle_mhz:

                # set frequency
                self.frequency_set(freq_mhz)

                # run sequence
                self.core_dma.playback_handle(handle_sequence)

                # add data to dataset
                with parallel:
                    self.update_dataset(freq_mhz, self.pmt_counter.fetch_count())
                    self.core.break_realtime()

        # reset devices
        self.dds_board.set_profile(0)
        self.dds_board.cfg_switches(0b1110)
        #fself.ttl_function_generator.off()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # tickle sequence
        with self.core_dma.record(_DMA_HANDLE_TICKLE_SEQUENCE):

            # start tickle
            with parallel:
                self.ttl_function_generator.on()
                delay_mu(self.time_fuction_generator_delay_mu)

            # read pmt
            self.pmt_gating_edge(self.time_tickle_us)

            # stop tickle
            #self.ttl_function_generator.off()


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up AOMs
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.core.break_realtime()

        # set up ttl for function generator trigger
        #self.ttl_function_generator.off()


    @rpc(flags={"async"})
    def frequency_set(self, freq_mhz):
        """
        Set the function generator to the desired frequency.
        """
        freq_set = self.fg.frequency(freq_mhz * 1e6)
        print('freq: {}'.format(freq_set))
        sleep(self.time_freq_delay_s)
        #print('freq set: {}'.format(freq_mhz))


    @rpc(flags={"async"})
    def update_dataset(self, freq_mhz, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('tickle_scan', [freq_mhz, pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        self.tickle_scan = np.array(self.tickle_scan)
        #self.fg.toggle(0)
        #self.ttl_function_generator.on()
        self.fg.deselect_device()
        #self.micromotion_compensation[:, 0] = float(self.micromotion_compensation[:, 0] / 2**16)
