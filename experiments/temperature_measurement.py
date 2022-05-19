import numpy as np
from artiq.experiment import *


class TemperatureMeasurement(EnvExperiment):
    """
    Measures ion fluorescence for a single detuning.
    """

    kernel_invariants = {
        "time_reset_mu",
        "dataset_name"
    }

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")             # always needed
        # experiment arguments
        self.setattr_argument("repetitions", NumberValue(default=1000, ndecimals=0, step=1, min=1, max=1000))
        # timing arguments
        self.setattr_argument("time_delay_us", NumberValue(ndecimals=2, step=1, min=1, max=1000))
        # 397nm
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")      # 397nm pump
        self.setattr_device("urukul0_ch1")      # 397nm probe
        self.setattr_argument("freq_probe_mhz", NumberValue(default=0, ndecimals=2, step=1, min=-100, max=100))
        self.setattr_argument("freq_pump_mhz", NumberValue(default=80, ndecimals=2, step=1, min=10, max=200))
        self.setattr_argument("time_pump_us", NumberValue(default=150, ndecimals=150, step=1, min=1, max=1000))
        self.setattr_argument("time_probe_us", NumberValue(default=50, ndecimals=50, step=1, min=1, max=1000))
        # PMT
        self.setattr_device("ttl0")             # PMT signal
        self.setattr_device("ttl4")             # PMT power
        self.setattr_device("ttl_counter0")     # PMT signal
        self.edge_method = "rising"             # read rising events
        # photodiode
        self.setattr_device("sampler0")
        self.photodiode_channel = 2

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # DDS devices
        self.dds_board = self.get_device('urukul0_cpld')
        self.dds_probe = self.get_device('urukul0_ch0')
        self.dds_pump = self.get_device('urukul0_ch1')
        # PMT devices
        self.ttl_signal = self.get_device('ttl0')
        self.ttl_power = self.get_device('ttl1')
        self.pmt_counter = self.get_device('ttl_counter0')
        # get edge gating method
        self.gate_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.edge_method))

        # set up datasets
        self.set_dataset("pmt_dataset", np.zeros(self.repetitions), broadcast=True)
        self.set_dataset("photodiode_reading", np.zeros(self.repetitions), broadcast=True)

        # convert time values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_mu * us)

        # convert dds values to machine units
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * 1e6)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)
        self.freq_probe_ftw = self.dds_probe.frequency_to_ftw(self.freq_probe_mhz * 1e6)
        self.ampl_probe_asf = self.dds_probe.amplitude_to_asf(0.5)
        # todo: set parameters

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
        handle = self.core_dma.get_handle("tmp")
        self.core.break_realtime()
        # read photodiode
        self.readPhotodiode()
        self.core.break_realtime()
        # run the experiment
        for trial_num in range(self.repetitions):
            # run pulse sequence from core DMA
            self.core_dma.playback_handle(handle)
            # record pmt counts to dataset
            self.mutate_dataset("pmt_dataset", trial_num, self.pmt_counter.fetch_counts())
            delay_mu(self.time_reset_mu)

    @kernel
    def DMArecord(self):
        """
        Record the 397nm sequence for a single data point onto core DMA.
        """
        with self.core_dma.record("tmp"):
            # pump on, probe off
            with parallel:
                self.dds_pump.cfg_sw(1)
                self.dds_probe.cfg_sw(0)
                delay_mu(self.time_pump_us)
            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_pump.cfg_sw(0)
                self.dds_probe.cfg_sw(1)
                self.pmt_counter.gate_edge(self.time_probe_us)
            with parallel:
                self.dds_pump.cfg_sw(0)
                self.dds_probe.cfg_sw(0)


    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()
        # set signal TTL to correct direction for PMT
        self.ttl_signal.input()
        # turn PMT on and wait for turn-on time
        self.ttl_power.on()
        delay_mu(560000)
        # initialize dds board
        self.dds_board.init()
        self.core.break_realtie()
        # initialize pump beam and set waveform
        self.dds_pump.init()
        self.core.break_realtime()
        self.dds_pump.set_mu(self.freq_pump_ftw, asf=self.ampl_pump_asf)
        self.dds_pump.set_att(10)
        self.core.break_realtime()
        # initialize probe beam and set waveform
        self.dds_probe.init()
        self.core.break_realtime()
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf)
        self.dds_probe.set_att(10)
        self.core.break_realtime()
        # set up sampler
        self.sampler0.init()
        self.core.break_realtime()
        self.sampler0.set_gain_mu(self.photodiode_channel, 1)

    @kernel
    def readPhotodiode(self):
        photodiode_buffer = np.zeros(self.photodiode_channel + 1)
        self.sampler0.sample(photodiode_buffer)
        self.core.break_realtime()
        photodiode_value = photodiode_buffer[self.photodiode_channel]
        self.mutate_dataset("photodiode_reading", 0, photodiode_value)
        self.core.break_realtime()

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # todo: upload data to labrad
        import labrad
        cxn = labrad.connect()
        print(cxn)
        pass
