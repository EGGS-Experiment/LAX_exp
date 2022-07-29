import numpy as np
from artiq.experiment import *

_DMA_HANDLE = "laser_scan"
# todo: make repump status more general
# todo: upload data to labrad
# todo: check synchronization of cycle with now_mu()
# todo: check scannable works correctly


class LaserScan(EnvExperiment):
    """
    729nm Laser Scan
    Measures ion fluorescence vs AOM detuning
    """

    #kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",            NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_delay_us",          NumberValue(default=100, ndecimals=2, step=1, min=1, max=1000))
        self.setattr_argument("time_probe_us",          NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000))

        # AOM DDS channels
        self.setattr_argument("dds_board_num",          NumberValue(default=0, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_qubit_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",    Scannable(default=RangeScan(65, 85, 5),
                                                                  global_min=60, global_max=110, global_step=1,
                                                                  unit="MHz", scale=1, ndecimals=1))

        # PMT
        self.setattr_argument("pmt_input_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",        EnumerationValue(["rising", "falling", "both"], default="rising"))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # DDS devices
        self.dds_board = self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))

        # convert dds values to machine units - probe
        self.ftw_to_frequency = 1e9 / (2**32 - 1)
        self.freq_probe_scan_mhz2 = list(self.freq_probe_scan_mhz)

        # convert dds values to machine units - everything else
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * MHz)
        self.freq_repump_ftw = self.dds_probe.frequency_to_ftw(self.freq_repump_mhz * MHz)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)

        # set up datasets
        self.set_dataset("temperature_measurement", [], broadcast=True)

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()
        self.core.break_realtime()

        # program pulse sequence onto core DMA
        self.DMArecord()
        self.core.break_realtime()

        # retrieve pulse sequence handle
        handle = self.core_dma.get_handle(_DMA_HANDLE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for freq_mhz in self.freq_probe_scan_mhz2:
            self.core.break_realtime()

            # set freq and ampl for probe
            freq_mu = self.dds_probe.frequency_to_ftw(freq_mhz)
            self.dds_probe.set_mu(freq_mu, asf=self.ampl_probe_asf)
            self.core.break_realtime()

            # run the experiment
            self.dds_repump.cfg_sw(1)
            for trial_num in range(self.repetitions):
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle)

                # update dataset
                with parallel:
                    self.update_dataset(freq_mhz, self.pmt_counter.fetch_count())
                    self.core.break_realtime()
                #delay_mu(self.time_delay_mu)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE):
            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_pump_states_on)
                delay_mu(self.time_pump_mu)
            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_probe_states_on)
                self.pmt_gating_edge(self.time_probe_mu)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # initialize dds board
        #self.dds_board.init()
        self.core.break_realtime()

        # initialize pump beam and set waveform
        #self.dds_pump.init()
        self.core.break_realtime()
        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf)
        self.core.break_realtime()
        self.dds_qubit.set_att(8 * dB)
        self.core.break_realtime()
        self.dds_qubit.cfg_sw(1)

    @rpc(flags={"async"})
    def update_dataset(self, freq_mhz, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('temperature_measurement', [freq_mhz, pmt_counts])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # data = self.get_dataset("temperature_measurement_repump_on")
        # pmt_data = data[2500:, 1]
        # pd_data = data[2500:, 2]
        #
        # print("repump on:")
        # print("\tpmt counts: {:.2f} +/- {:.2f}".format(np.mean(pmt_data), np.std(pmt_data)))
        # print("\tphotodiode value: {:.2f} +/- {:.2f}".format(np.mean(pd_data), np.std(pd_data)))
