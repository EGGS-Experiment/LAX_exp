import numpy as np
from artiq.experiment import *

#_DMA_HANDLE = "temperature_measurement"
_DMA_HANDLE_ON = "temperature_measurement_on"
_DMA_HANDLE_OFF = "temperature_measurement_off"
# todo: synchronize photodiode with pmt
# todo: frequency 75 +/- 10mhz, 5 points total
# todo: 40mV from photodiode
# todo: make repump status more general
# todo: upload data to labrad
# todo: set attenuations as arguments
# todo: convert reading to actual volts, will be 400mV
# todo: use flags to speed things up
# todo: reset all dds channels to original state
# todo: check synchronization of cycle
# todo: implement scan


class TemperatureMeasurement(EnvExperiment):
    """
    Temperature Measurement
    Measures ion fluorescence for a single detuning.
    """

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

        # AOM DDSs
        self.setattr_argument("dds_board_num",          NumberValue(default=1, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_probe_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_pump_channel",       NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_repump_channel",     NumberValue(default=2, ndecimals=0, step=1, min=0, max=3))

        # AOM parameters
        self.setattr_argument("freq_probe_mhz",         NumberValue(default=70, ndecimals=3, step=1, min=-100, max=100))
        self.setattr_argument("freq_pump_mhz",          NumberValue(default=140, ndecimals=3, step=1, min=10, max=200))
        self.setattr_argument("freq_repump_mhz",        NumberValue(default=110, ndecimals=3, step=1, min=10, max=200))

        # PMT
        self.setattr_argument("pmt_input_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        #self.setattr_argument("pmt_power_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",        EnumerationValue(["rising", "falling", "both"], default="rising"))

        # photodiode
        # todo: see if we can move sampler setattr to prepare
        self.setattr_device("sampler0")
        self.setattr_argument("photodiode_channel",     NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("photodiode_gain",        NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # DDS devices
        self.dds_board = self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_channel))

        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))
        #self.ttl_signal = self.get_device("ttl{:d}".format(self.pmt_input_channel))
        #self.ttl_power = self.get_device("ttl{:d}".format(self.pmt_power_channel))

        # convert time values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # convert dds values to machine units
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * 1e6)
        self.freq_probe_ftw = self.dds_probe.frequency_to_ftw(self.freq_probe_mhz * 1e6)
        self.freq_repump_ftw = self.dds_probe.frequency_to_ftw(self.freq_repump_mhz * 1e6)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)
        self.ampl_probe_asf = self.dds_probe.amplitude_to_asf(0.5)
        self.ampl_repump_asf = self.dds_probe.amplitude_to_asf(0.5)

        # frequency scan range
        self.freq_probe_range = [self.dds_pump.frequency_to_ftw(freq_tmp * 1e6) for freq_tmp in [75]]

        # get DDS board switch states
        # todo: sort out 866 repump generally, as well as 854 always on
        # self.dds_switch_pump_states = 0b0000 | (0b1 << self.dds_pump_channel)
        # self.dds_switch_probe_states = 0b0000 | (0b1 << self.dds_probe_channel)
        self.dds_switch_pump_states_on = 0b1100 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_on = 0b1100 | (0b1 << self.dds_probe_channel)
        self.dds_switch_pump_states_off = 0b1000 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_off = 0b1000 | (0b1 << self.dds_probe_channel)

        # set up datasets
        # todo: do a better way
        #self.set_dataset("temperature_measurement", np.zeros([len(self.freq_probe_range), self.repetitions * 2, 4]), broadcast=True)
        self.set_dataset("temperature_measurement", np.zeros([len(self.freq_probe_range) * self.repetitions * 2, 4]), broadcast=True)

        # adc
        self.adc = self.sampler0
        self.adc_mu_to_volts = (10 ** (1 - self.photodiode_gain)) / (2 ** 15)

    @kernel
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

        # set up photodiode buffer
        photodiode_buffer = [0] * 8

        # retrieve pulse sequence handle
        #handle = self.core_dma.get_handle(_DMA_HANDLE)
        handle_on = self.core_dma.get_handle(_DMA_HANDLE_ON)
        handle_off = self.core_dma.get_handle(_DMA_HANDLE_OFF)
        self.core.break_realtime()

        # run the experiment (repump on)
        for i in range(len(self.freq_probe_range)):

            # set ftw
            self.dds_probe.set_mu(self.freq_probe_range[i], asf=self.ampl_probe_asf)
            self.core.break_realtime()

            for trial_num in range(self.repetitions):
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_on)

                # record data
                # maybe: try parallel?
                self.adc.sample_mu(photodiode_buffer)
                #self.core.break_realtime()
                self.mutate_dataset("temperature_measurement", trial_num,
                                    [self.freq_probe_range[i],
                                     1,
                                     self.pmt_counter.fetch_count(),
                                     photodiode_buffer[self.photodiode_channel]
                                     ])
                # self.mutate_dataset("temperature_measurement", ((i, i + 1), (trial_num, trial_num + 1)),
                #                     [self.freq_probe_range[i],
                #                      1,
                #                      self.pmt_counter.fetch_count(),
                #                      photodiode_buffer[self.photodiode_channel]
                #                      ])
                self.core.break_realtime()
                #delay_mu(self.time_delay_mu)

            # run the experiment (repump off)
            self.dds_repump.cfg_sw(0)
            for trial_num in range(self.repetitions, 2 * self.repetitions):
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_off)

                # record data
                # maybe: try parallel?
                self.adc.sample_mu(photodiode_buffer)
                #self.core.break_realtime()
                self.mutate_dataset("temperature_measurement", trial_num,
                                    [self.freq_probe_range[i],
                                     0,
                                     self.pmt_counter.fetch_count(),
                                     photodiode_buffer[self.photodiode_channel]
                                     ])
                # self.mutate_dataset("temperature_measurement", ((i, i + 1), (trial_num, trial_num + 1)),
                #                     [self.freq_probe_range[i],
                #                      0,
                #                      self.pmt_counter.fetch_count(),
                #                      photodiode_buffer[self.photodiode_channel]
                #                      ])
                #delay_mu(self.time_delay_mu)
                self.core.break_realtime()

            self.dds_repump.cfg_sw(1)
            self.dds_pump.cfg_sw(1)

    @kernel
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

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set signal TTL to correct direction for PMT
        self.core.break_realtime()
        #self.ttl_signal.input()

        # turn PMT on and wait for turn-on time
        #self.ttl_power.on()
        #delay_mu(560000)

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