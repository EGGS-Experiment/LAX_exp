import numpy as np
from artiq.experiment import *

_DMA_HANDLE = "temperature_measurement"


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
        self.setattr_argument("repetitions",            NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))

        # timing
        self.setattr_argument("time_delay_us",          NumberValue(default=100, ndecimals=2, step=1, min=1, max=1000))

        # AOM DDSs
        self.setattr_argument("dds_board_num",          NumberValue(default=1, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_probe_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("dds_pump_channel",       NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))

        # AOM parameters
        self.setattr_argument("freq_probe_mhz",         NumberValue(default=110, ndecimals=2, step=1, min=-100, max=100))
        self.setattr_argument("time_probe_us",          NumberValue(default=500, ndecimals=50, step=1, min=1, max=1000))
        self.setattr_argument("freq_pump_mhz",          NumberValue(default=140, ndecimals=2, step=1, min=10, max=200))
        self.setattr_argument("time_pump_us",           NumberValue(default=500, ndecimals=150, step=1, min=1, max=1000))

        # PMT
        self.setattr_argument("pmt_input_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_power_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",        EnumerationValue(["rising", "falling", "both"], default="rising"))

        # photodiode
        self.setattr_device("sampler0")
        self.setattr_argument("photodiode_channel",     NumberValue(default=2, ndecimals=0, step=1, min=0, max=7))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # DDS devices
        self.dds_board = self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))

        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))
        self.ttl_signal = self.get_device("ttl{:d}".format(self.pmt_input_channel))
        self.ttl_power = self.get_device("ttl{:d}".format(self.pmt_power_channel))

        # set up datasets
        self.set_dataset("pmt_dataset", np.zeros(self.repetitions), broadcast=True)
        self.set_dataset("photodiode_reading", 0, broadcast=True)

        # convert time values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_us * us)

        # convert dds values to machine units
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * 1e6)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)
        self.freq_probe_ftw = self.dds_probe.frequency_to_ftw(self.freq_probe_mhz * 1e6)
        self.ampl_probe_asf = self.dds_probe.amplitude_to_asf(0.5)

        # get DDS board switch states
        self.dds_switch_pump_states = 0b0000 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states = 0b0000 | (0b1 << self.dds_probe_channel)

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
        # retrieve pulse sequence handle
        handle = self.core_dma.get_handle(_DMA_HANDLE)
        self.core.break_realtime()
        # read photodiode
        self.readPhotodiode()
        self.core.break_realtime()
        # run the experiment
        for trial_num in range(self.repetitions):
            # run pulse sequence from core DMA
            self.core_dma.playback_handle(handle)
            # record pmt counts to dataset
            self.mutate_dataset("pmt_dataset", trial_num, self.pmt_counter.fetch_count())
            delay_mu(self.time_delay_mu)

    @kernel
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE):
            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_pump_states)
                delay_mu(self.time_pump_mu)
            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(self.dds_switch_probe_states)
                self.pmt_gating_edge(self.time_probe_mu)

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set signal TTL to correct direction for PMT
        self.core.break_realtime()
        self.ttl_signal.input()

        # turn PMT on and wait for turn-on time
        self.ttl_power.on()
        delay_mu(560000)

        # initialize dds board
        self.dds_board.init()
        self.core.break_realtime()

        # initialize pump beam and set waveform
        self.dds_pump.init()
        self.core.break_realtime()
        self.dds_pump.set_mu(self.freq_pump_ftw, asf=self.ampl_pump_asf)
        self.dds_pump.set_att(12 * dB)
        self.core.break_realtime()
        self.dds_pump.cfg_sw(1)

        # initialize probe beam and set waveform
        self.dds_probe.init()
        self.core.break_realtime()
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf)
        self.dds_probe.set_att(3 * dB)
        self.core.break_realtime()

        # set up sampler
        self.sampler0.init()
        self.core.break_realtime()
        self.sampler0.set_gain_mu(self.photodiode_channel, 1)

    @kernel
    def readPhotodiode(self):
        # get sampler value
        photodiode_buffer = [0.0] * 8
        self.sampler0.sample(photodiode_buffer)
        self.core.break_realtime()

        # record value
        photodiode_value = photodiode_buffer[self.photodiode_channel]
        self.set_dataset("photodiode_reading", photodiode_value, broadcast=True)
        self.core.break_realtime()

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        photodiode_value = self.get_dataset("photodiode_reading")
        # print('photodiode value:')
        # print(photodiode_value)

        pmt_counts = self.get_dataset("pmt_dataset")
        print("pmt counts:")
        print(pmt_counts)

        # todo: upload data to labrad
        # import labrad
        # cxn = labrad.connect()
        # print(cxn)
        # pass
