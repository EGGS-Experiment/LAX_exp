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
        self.setattr_argument("freq_detuning_mhz", NumberValue(default=0, ndecimals=2, step=1, min=-100, max=100))
        self.setattr_argument("time_pump_us", NumberValue(default=150, ndecimals=150, step=1, min=1, max=1000))
        self.setattr_argument("time_probe_us", NumberValue(default=50, ndecimals=50, step=1, min=1, max=1000))
        # PMT
        self.setattr_device("ttl0")             # PMT signal
        self.setattr_device("ttl4")             # PMT power

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # todo: re-alias devices to have more sensible names
        # get TTLs
        self.ttl_signal = self.get_device('ttl{:d}'.format(self.signal_ttl))
        self.ttl_power = self.get_device('ttl{:d}'.format(self.power_ttl))

        # check TTLs are input
        if (self.ttl_signal.__class__.__name__ != 'TTLInOut') or (self.ttl_trigger.__class__.__name__ != 'TTLInOut'):
            raise Exception('Error: TTLs must be input.')
        # check linetrigger and PMT are on different TTLs
        elif self.trigger_status and (self.ttl_signal == self.ttl_trigger):
            raise Exception('Error: Trigger and PMT TTLs cannot be the same.')
        # check edge gating input valid
        elif self.edge_method not in ('rising', 'falling', 'both'):
            raise Exception('Error: invalid edge method.')
        # get edge gating method
        self.gate_edge = getattr(self.ttl_signal, 'gate_{:s}_mu'.format(self.edge_method))

        # set up dataset
        #date = datetime.now()
        #self.dataset_name = 'pmt_{:s}_{:02d}_{:02d}_{:02d}:{:02d}'.format(str(date.year), date.month, date.day, date.hour, date.minute)
        self.dataset_name = "tmpexpres"
        self.set_dataset(self.dataset_name, np.zeros(self.num_bins), broadcast=True)

        # convert any values to machine units
        self.time_pump_mu = self.core.seconds_to_mu(self.time_pump_us * us)
        self.time_probe_mu = self.core.seconds_to_mu(self.time_probe_us * us)
        self.time_delay_mu = self.core.seconds_to_mu(self.time_delay_mu * us)
        # todo
        self.freq_detuning_mhz = None

    @kernel
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()
        # program the laser sequence
        self.core.DMArecord()
        self.core.break_realtime()
        # todo: get the handle for the laser sequence
        # prepare the PMT
        self.PMTprepare()
        self.core.break_realtime()
        # run the experiment
        for trial_num in range(self.repetitions):
            self.laserPrepare()
            self.science(trial_num)
            self.core.break_realtime()

    @kernel
    def DMArecord(self):
        """
        Record the 397nm sequence for a single data point
        onto core DMA.
        """
        with self.core_dma.record("tmp"):
            pass
            # todo: set the probe beam switch
            # todo:
        # todo: see if we can set handle here

    @kernel
    def PMTprepare(self):
        """
        Prepare the PMT for the experiment.
        """
        # set signal TTL to correct direction
        self.ttl_signal.input()
        # switch the PMT on
        self.ttl_power.on()

    @kernel
    def laserPrepare(self):
        """
        Set the lasers to the correct frequencies.
        """
        # todo: turn all lasers off
        # todo: set the probe beam waveform
        # todo: set the probe beam attenuation
        # todo: set the pump beam waveform
        # todo: set the pump beam attenuation
        pass

    @kernel
    def science(self, trial_num):
        """
        Actually run the sequence and do science.
        """
        # todo: call dma sequence
        # todo: record pmt counts
        # todo: switch off pmt
        pass


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # todo: try to upload to labrad
        pass
