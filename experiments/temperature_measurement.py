from artiq.experiment import *
import numpy as np


class TemperatureMeasurement(EnvExperiment):
    """
    Measures ion fluorescence for a single detuning.
    """

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
        self.edge_method = "rising"             # read rising events

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

        # check TTLs are input
        if self.ttl_signal.__class__.__name__ != 'TTLInOut':
            raise Exception('Error: TTL must be input.')
        # check edge gating input valid
        elif self.edge_method not in ('rising', 'falling', 'both'):
            raise Exception('Error: invalid edge method.')
        # get edge gating method
        self.gate_edge = getattr(self.ttl_signal, 'gate_{:s}_mu'.format(self.edge_method))

        # set up dataset
        #date = datetime.now()
        #self.dataset_name = 'pmt_{:s}_{:02d}_{:02d}_{:02d}:{:02d}'.format(str(date.year), date.month, date.day, date.hour, date.minute)
        self.dataset_name = "tmp_exp_res"
        self.set_dataset(self.dataset_name, np.zeros(self.num_bins), broadcast=True)

        # convert time values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_mu * us)

        # convert dds values to machine units
        self.freq_pump_ftw = self.dds_pump.frequency_to_ftw(self.freq_pump_mhz * 1e6)
        self.ampl_pump_asf = self.dds_pump.amplitude_to_asf(0.5)
        self.freq_probe_ftw = self.dds_probe.frequency_to_ftw(self.freq_probe_mhz * 1e6)
        self.ampl_probe_asf = self.dds_probe.amplitude_to_asf(0.5)

    @kernel
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()
        # prepare devices
        self.prepareDevices()
        self.core.break_realtime()
        # run the experiment
        for trial_num in range(self.repetitions):
            # turn on pump beam and turn off probe beam (if it wasn't already off)
            with parallel:
                self.dds_pump.cfg_sw(1)
                self.dds_probe.cfg_sw(0)
            # wait to cool
            delay_mu(self.time_pump_us)
            # turn on probe beam and turn off pump beam
            with parallel:
                self.dds_pump.cfg_sw(0)
                self.dds_probe.cfg_sw(1)
            # record pmt counts
            self.mutate_dataset(self.dataset_name, i, self.ttl_signal.count(self.gate_edge(self.time_probe_mu)))
            # turn off probe beam
            self.dds_probe.cfg_sw(0)
            # turn on pump
            self.dds_pump.cfg_sw(1)
            delay_mu(self.time_reset_mu)

    @kernel
    def DMArecord(self):
        """
        Record the 397nm sequence for a single data point
        onto core DMA.
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
                self.dds_probe.cfg_sw(0)
            # todo: implement edge counter to read out

    @kernel
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()
        # set signal TTL to correct direction for PMT
        self.ttl_signal.input()
        # set pump beam waveform
        self.dds_pump.set_mu(self.freq_pump_ftw, asf=self.ampl_pump_asf)
        self.core.break_realtime()
        # set probe beam waveform
        self.dds_probe.set_mu(self.freq_probe_ftw, asf=self.ampl_probe_asf)
        self.core.break_realtime()

    @kernel
    def prepareTrial(self):
        """
        Prepare an individual trial.
        Set the lasers to their correct frequencies and turn on the PMT
        """
        pass

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # todo: try to upload to labrad
        pass
