import numpy as np
from artiq.experiment import *
from artiq.language import *

#_DMA_HANDLE = "temperature_measurement"
_DMA_HANDLE_ON = "temperature_measurement_on"
_DMA_HANDLE_OFF = "temperature_measurement_off"
# todo: make repump status more general
# todo: upload data to labrad
# todo: check synchronization of cycle with now_mu()
# todo: check scannable works correctly


class TemperatureMeasurement(EnvExperiment):
    """
    Temperature Measurement
    Measures ion fluorescence for a single detuning.
    """

    #kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")             # always needed
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",            NumberValue(default=7500, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_delay_us",          NumberValue(default=100, ndecimals=2, step=1, min=1, max=1000))
        self.setattr_argument("time_probe_us",          NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000))
        self.setattr_argument("time_pump_us",           NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000))

        # AOM DDS channels
        self.setattr_argument("dds_board_num",          NumberValue(default=1, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_probe_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_pump_channel",       NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_repump_channel",     NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",    Scannable(default=RangeScan(65, 85, 5),
                                                                  global_min=60, global_max=110, global_step=1,
                                                                  unit = "MHz", scale = 1, ndecimals = 1))

        # AOM parameters
        self.setattr_argument("freq_pump_mhz",          NumberValue(default=140, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_repump_mhz",        NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))

        # PMT
        self.setattr_argument("pmt_input_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",        EnumerationValue(["rising", "falling", "both"], default="rising"))

        # photodiode
        self.setattr_argument("photodiode_channel",     NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("photodiode_gain",        NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # DDS devices
        self.dds_board = self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_channel))

        # convert dds values to machine units
        self.ftw_to_frequency = 1e9 / (2**32 - 1)
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * MHz)
        self.freq_repump_ftw = self.dds_probe.frequency_to_ftw(self.freq_repump_mhz * MHz)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)
        self.ampl_probe_asf = self.dds_probe.amplitude_to_asf(0.5)
        self.ampl_repump_asf = self.dds_probe.amplitude_to_asf(0.5)

        # get DDS board switch states, on/off signifies repump on/off
        self.dds_switch_pump_states_on = 0b1100 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_on = 0b1100 | (0b1 << self.dds_probe_channel)
        self.dds_switch_pump_states_off = 0b1000 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_off = 0b1000 | (0b1 << self.dds_probe_channel)

        # ADC
        self.adc = self.get_device("sampler0")
        self.adc_buffer = np.zeros(8, dtype=int)
        self.adc_mu_to_volts = (10 ** (1 - self.photodiode_gain)) / (2 ** 15)

        # set up datasets
        self.set_dataset("temperature_measurement_repump_on", np.zeros([len(self.freq_probe_scan_mhz) * self.repetitions, 3]), broadcast=True)
        self.set_dataset("temperature_measurement_repump_off", np.zeros([len(self.freq_probe_scan_mhz) * self.repetitions, 3]), broadcast=True)

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        #self.prepareDevices()
        self.core.break_realtime()

        # program pulse sequence onto core DMA
        self.DMArecord()
        self.core.break_realtime()

        # retrieve pulse sequence handle
        handle_on = self.core_dma.get_handle(_DMA_HANDLE_ON)
        handle_off = self.core_dma.get_handle(_DMA_HANDLE_OFF)
        self.core.break_realtime()


        # MAIN SEQUENCE
        for freq_mhz in self.freq_probe_scan_mhz:
            set freq and ampl for probe
            freq_mu = self.dds_probe.frequency_to_ftw(freq_mhz)
            #self.dds_probe.set_mu(freq_mu, asf=self.ampl_probe_asf)
            self.core.break_realtime()

            # run the experiment (repump on)
            self.dds_repump.cfg_sw(1)
            for trial_num in range(self.repetitions):
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_on)
                #self.core.break_realtime()

                # get data
                with parallel:
                    count = self.pmt_counter.fetch_count()
                    self.adc.sample_mu(self.adc_buffer)

                # update dataset
                self.update_dataset(trial_num, freq_mhz, 1, count, self.adc_buffer[self.photodiode_channel])
                self.core.break_realtime()
                #delay_mu(self.time_delay_mu)

            # run the experiment (repump off)
            self.dds_repump.cfg_sw(0)
            for trial_num in range(self.repetitions):
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_off)
                #self.core.break_realtime()

                # record data
                with parallel:
                    count = self.pmt_counter.fetch_count()
                    self.adc.sample_mu(self.adc_buffer)

                # update dataset
                self.update_dataset(trial_num, freq_mhz, 0, count, self.adc_buffer[self.photodiode_channel])
                self.core.break_realtime()
                #delay_mu(self.time_delay_mu)

            # after sequence, set all dds channels to trapping state
            self.dds_repump.cfg_sw(1)
            self.dds_pump.cfg_sw(1)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # tmp remove: bring back parallel
        with self.core_dma.record(_DMA_HANDLE_ON):
            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_pump_states_on)
                delay_mu(self.time_pump_mu)
            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_probe_states_on)
                self.pmt_gating_edge(self.time_probe_mu)

        self.core.break_realtime()
        with self.core_dma.record(_DMA_HANDLE_OFF):
            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_pump_states_off)
                delay_mu(self.time_pump_mu)
            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_probe_states_off)
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
        self.dds_pump.set_mu(self.freq_pump_ftw, asf=self.ampl_pump_asf)
        self.core.break_realtime()
        self.dds_pump.set_att(13 * dB)
        self.core.break_realtime()
        self.dds_pump.cfg_sw(1)

        # initialize probe beam and set waveform
        #self.dds_probe.init()
        self.core.break_realtime()
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf)
        self.dds_probe.set_att(8 * dB)
        self.core.break_realtime()
        self.dds_pump.cfg_sw(0)

        # initialize repump beam and set waveform
        #self.dds_repump.init()
        self.core.break_realtime()
        self.dds_repump.set_mu(self.freq_repump_ftw, asf=self.ampl_repump_asf)
        self.dds_repump.set_att(6.5 * dB)
        self.core.break_realtime()
        self.dds_repump.cfg_sw(1)

        # set up sampler
        #self.adc.init()
        self.core.break_realtime()
        self.adc.set_gain_mu(self.photodiode_channel, self.photodiode_gain)

    @rpc(flags={"async"})
    def thkim(self, yz):
        """
        Records values via rpc to minimize kernel overhead.
        """
        print(yz)

    @rpc(flags={"async"})
    def update_dataset(self, i, freq_mu, repump_status, pmt_counts, sampler_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        if repump_status == 1:
            self.mutate_dataset("temperature_measurement_repump_on", i, [freq_mu * self.ftw_to_frequency, pmt_counts, sampler_mu * self.adc_mu_to_volts])
        else:
            self.mutate_dataset("temperature_measurement_repump_off", i, [freq_mu * self.ftw_to_frequency, pmt_counts, sampler_mu * self.adc_mu_to_volts])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        data = self.get_dataset("temperature_measurement")
        pmt_data = data[2500:, 2]
        pd_data = data[2500:, 3]

        print("repump off:")
        print("\tpmt counts: {:.2f} +/- {:.2f}".format(np.mean(pmt_data), np.std(pmt_data)))
        print("\tphotodiode value: {:.2f} +/- {:.2f}".format(np.mean(pd_data), np.std(pd_data)))
