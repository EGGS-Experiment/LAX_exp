import numpy as np
from artiq.experiment import *

_DMA_HANDLE_ON = "temperature_measurement_on"
_DMA_HANDLE_OFF = "temperature_measurement_off"
# todo: make repump status more general
# todo: upload data to labrad
# todo: check synchronization of cycle with now_mu()
# todo: check scannable works correctly
# todo: set up ion calibration properly
# todo: fix ADC
# todo: set attenuations problem
# todo: remove freq_probe_scan_mhz2 and fix scannables problem


class TemperatureMeasurement(EnvExperiment):
    """
    Temperature Measurement - Transpose
    Measures ion fluorescence for a single detuning.
    Sweeps frequency within each trial.
    """

    #kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")             # always needed
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",            NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_probe_us",          NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",    Scannable(default=RangeScan(70, 146, 20, randomize=True),
                                                                  global_min=80, global_max=140, global_step=1,
                                                                  unit="MHz", scale=1, ndecimals=1))

        # photodiode
        self.setattr_argument("photodiode_channel",     NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("photodiode_gain",        NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =              self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_pump_mu =                 self.core.seconds_to_mu(self.time_doppler_cooling * us)
        self.time_probe_mu =                self.core.seconds_to_mu(self.time_probe_us * us)

        # DDS devices
        self.dds_board =                    self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe =                    self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                     self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump =                   self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_channel))

        # convert dds values to machine units - probe
        self.ftw_to_frequency =             1e9 / (2**32 - 1)
        self.freq_probe_scan_mhz2 =         list(self.freq_probe_scan_mhz)
        self.freq_probe_scan_ftw =          [self.dds_probe.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz2]

        # convert dds values to machine units - everything else
        self.freq_repump_ftw =              self.dds_probe.frequency_to_ftw(self.freq_repump_mhz * MHz)
        self.ampl_pump_asf =                self.dds_pump.amplitude_to_asf(0.5)
        self.ampl_probe_asf =               self.dds_probe.amplitude_to_asf(0.5)
        self.ampl_repump_asf =              self.dds_probe.amplitude_to_asf(0.5)

        # get DDS board switch states, on/off signifies repump on/off
        self.dds_switch_pump_states_on =    0b0100 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_on =   0b0100 | (0b1 << self.dds_probe_channel)
        self.dds_switch_pump_states_off =   0b0000 | (0b1 << self.dds_pump_channel)
        self.dds_switch_probe_states_off =  0b0000 | (0b1 << self.dds_probe_channel)

        # ADC
        self.adc = self.get_device("sampler0")
        self.adc_buffer = np.zeros(8, dtype=int)
        self.adc_mu_to_volts = (10 ** (1 - self.photodiode_gain)) / (2 ** 15)

        # set up datasets
        self.set_dataset("temperature_measurement", [])
        self.setattr_dataset("temperature_measurement")
        self.set_dataset("temperature_measurement_processed", np.zeros([len(self.freq_probe_scan_mhz), 3]))
        self.setattr_dataset("temperature_measurement_processed")

        # attenuations:
        # original
        #self.att_probe = [6.5, 8.5, 10, 11.5, 12.5, 13, 13, 13, 12.5, 11, 8.5, 6]
        # below is for 30uW @729nm setting
        #self.att_probe = [25.0, 24.5, 24.0, 22.0, 23.5, 24.5, 24.5, 24.0, 24.0, 23.5, 22.5]
        # below is for 30uW @729nm setting, but with more points (4MHz steps, 70 to 146 MHz)
        #self.att_probe = [24.5, 24.0, 24.0, 24.0, 24.5, 25.0, 24.5, 24.0, 22.0, 23.5, 24.5, 24.5, 24.0, 24.0, 23.5, 22.5, 21.5, 19.5, 16.5, 13.0]
        # below is for 5uW @729nm setting
        att_freqs = np.linspace(70, 146, 20)
        att_vals = np.array([27.0, 26.5, 26.5, 26.5, 27.0, 27.0, 26.5, 26.0, 24.0, 25.5, 26.5, 27.0, 26.5, 26.0, 25.5, 24.5, 23.5, 22.0, 19.0, 15.5])
        att_dict = dict(np.concatenate([[att_freqs], [att_vals]]).transpose())
        self.att_probe = [att_dict[freq] for freq in self.freq_probe_scan_mhz2]
        self.att_probe = [np.int32(0xFF) - np.int32(round(att_dB * 8)) for att_dB in self.att_probe]
        self.att_reg = 0x00000000

        # tmp remove
        self.setattr_device('urukul1_ch3')

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
        handle_on = self.core_dma.get_handle(_DMA_HANDLE_ON)
        handle_off = self.core_dma.get_handle(_DMA_HANDLE_OFF)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):
            for i in range(len(self.freq_probe_scan_mhz2)):
                self.core.break_realtime()

                # set freq and ampl for probe
                freq_mhz = self.freq_probe_scan_mhz2[i]
                self.dds_probe.set_mu(self.freq_probe_scan_ftw[i], asf=self.ampl_probe_asf)
                self.core.break_realtime()

                # set att for probe
                att_reg_tmp = self.att_reg | (self.att_probe[i] << (self.dds_probe_channel * 8))
                self.dds_board.set_all_att_mu(att_reg_tmp)
                self.core.break_realtime()

                # run the trial (repump on)
                self.dds_repump.cfg_sw(1)
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_on)

                # update dataset
                with parallel:
                    self.update_dataset(freq_mhz, 1, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

                # run the experiment (repump off)
                self.dds_repump.cfg_sw(0)
                # run pulse sequence from core DMA
                self.core_dma.playback_handle(handle_off)
                with parallel:
                    self.update_dataset(freq_mhz, 0, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

            # after sequence, set all dds channels to trapping state
            self.dds_repump.cfg_sw(1)
            self.dds_pump.cfg_sw(1)
            self.dds_probe.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
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

        # get current attenuation register status
        self.att_reg = np.int32(self.dds_board.get_att_mu())
        self.att_reg &= ~(0xFF << (8 * self.dds_probe_channel))

        # initialize repump beam and set waveform
        self.core.break_realtime()
        self.dds_repump.set_mu(self.freq_repump_ftw, asf=self.ampl_repump_asf)
        self.dds_repump.cfg_sw(1)
        self.core.break_realtime()

        # set up sampler
        self.core.break_realtime()
        self.adc.set_gain_mu(self.photodiode_channel, self.photodiode_gain)

    @rpc(flags={"async"})
    def update_dataset(self, freq_mhz, repump_status, pmt_counts, sampler_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('temperature_measurement', [freq_mhz, repump_status, pmt_counts, sampler_mu * self.adc_mu_to_volts])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.temperature_measurement = np.array(self.temperature_measurement)

        # get sorted x-values (frequency)
        freq_list_mhz = sorted(set(self.temperature_measurement[:, 0]))

        # collate results
        collated_results = {
            # two lists, one for w/repump (index 1), and one for w/o repump (index 0)
            freq: [[], []]
            for freq in freq_list_mhz
        }
        for freq_mhz, repump_status, pmt_counts, sampler_volts in self.temperature_measurement:
            collated_results[freq_mhz][int(repump_status)].append(pmt_counts)

        # process counts for mean and std and put into processed dataset
        for i, (freq_mhz, count_list) in enumerate(collated_results.items()):
            self.temperature_measurement_processed[i] = np.array([freq_mhz, np.mean(count_list[1]) - np.mean(count_list[0]), np.std(count_list[1])])
