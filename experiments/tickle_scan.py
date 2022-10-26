import numpy as np
from time import sleep
from artiq.experiment import *

_DMA_HANDLE_TICKLE_SEQUENCE = "tickle_scan_sequence"

import labrad
from os import environ


class TickleScan(EnvExperiment):
    """
    Tickle Scan
    Scans the RF from a function generator and measures fluorescence.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",
        "ttl_channel_function_generator"

        "time_profileswitch_delay_us",
        "time_redist_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",

        "dds_board_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",

        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_repump_cooling_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
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
        self.setattr_argument("repetitions",                        NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))

        # tickle values
        self.setattr_argument("time_tickle_us",                     NumberValue(default=300, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("freq_tickle_mhz",                    Scannable(
                                                                        default=RangeScan(0.9, 1.2, 301),
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

        # convert time values to machine units
        self.time_profileswitch_delay_mu =                          self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_cooling_mu =                                      self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                                      self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_redist_mu =                                       self.core.seconds_to_mu(self.time_redist_us * us)

        # DDS devices
        self.dds_board =                                            self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                      self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                            self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                             self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                   self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))

        # convert frequency to ftw
        self.freq_redist_ftw =                                      self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                                self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =                              self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                      self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                                self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                              self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

        # tickle
        self.time_tickle_us =                                       self.core.seconds_to_mu(self.time_tickle_us * us)

        # prepare labrad devices
        self.fg.select_device()
        # todo: set up triggering

        # set up datasets
        self.set_dataset("tickle_scan", [])
        self.setattr_dataset("tickle_scan")


    #@kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        self.handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_TICKLE_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for freq_mhz in self.freq_tickle_mhz:

            # set frequency
            self.frequency_set(freq_mhz)

            # run experiment
            self.run_repetition(freq_mhz)

        # finish experiment
        self.run_reset()


    @kernel(flags={"fast-math"})
    def run_repetition(self, freq_mhz):
        """
        todo: document
        """
        # repeat experiment
        for trial_num in range(self.repetitions):
            # run sequence
            self.core_dma.playback_handle(self.handle_sequence)

            # add data to dataset
            with parallel:
                self.update_dataset(freq_mhz, self.pmt_counter.fetch_count())
                self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def run_reset(self):
        """
        Reset devices for normal operation.
        """
        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # cooling sequence
        with self.core_dma.record(_DMA_HANDLE_TICKLE_SEQUENCE):

            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # do cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # todo: check if spin depolarization needed
            # do spin depolarization using probe
            # self.dds_board.cfg_switches(0b0101)
            # delay_mu(self.time_redist_mu)
            # self.dds_board.cfg_switches(0b0100)

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # trigger function generator
            with parallel:
                self.ttl_function_generator.on()
                # todo: delay time it takes for function generator to turn on

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
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def frequency_set(self, freq_mhz):
        """
        Set the function generator to the desired frequency.
        """
        self.fg.frequency(freq_mhz * 1e6)
        print('freq set: {}'.format(freq_mhz))


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
        #self.micromotion_compensation[:, 0] = float(self.micromotion_compensation[:, 0] / 2**16)
        pass
