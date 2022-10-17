import numpy as np
from artiq.experiment import *

_DMA_HANDLE_ON = "temperature_measurement_on"
_DMA_HANDLE_OFF = "temperature_measurement_off"


class TemperatureMeasurement(EnvExperiment):
    """
    Temperature Measurement
    Measures ion fluorescence for a single detuning and sweeps frequency within each trial.
    """

    #kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_doppler_cooling_us",
        "time_readout_us",

        "dds_board_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
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
        self.setattr_argument("repetitions",            NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_probe_us",          NumberValue(default=100, ndecimals=5, step=1, min=1, max=1000))

        # probe parameters
        self.setattr_argument("ampl_probe_pct",         NumberValue(default=50, ndecimals=3, step=1, min=10, max=100))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",    Scannable(default=RangeScan(70, 146, 20, randomize=True),
                                                                  global_min=10, global_max=200, global_step=1,
                                                                  unit="MHz", scale=1, ndecimals=5))

        # photodiode
        self.setattr_argument("photodiode_channel",     NumberValue(default=0, ndecimals=0, step=1, min=0, max=7))
        self.setattr_argument("photodiode_gain",        NumberValue(default=1, ndecimals=0, step=1, min=0, max=3))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                              self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                  self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_probe_mu =                            self.core.seconds_to_mu(self.time_probe_us * us)

        # DDS devices
        self.dds_board =                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_probe =                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))

        # convert dds values to machine units - probe
        self.ftw_to_frequency =                         1e9 / (2**32 - 1)
        self.freq_probe_scan_ftw =                      [self.dds_probe.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz]

        # convert dds values to machine units - everything else
        self.freq_pump_cooling_ftw =                    self.dds_probe.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_repump_cooling_ftw =                  self.dds_probe.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)

        self.ampl_probe_asf =                           self.dds_probe.amplitude_to_asf(self.ampl_probe_pct / 100)
        self.ampl_pump_cooling_asf =                    self.dds_probe.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_repump_cooling_asf =                  self.dds_probe.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

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
        self.att_probe = [att_dict[freq] for freq in self.freq_probe_scan_mhz]
        self.att_probe = [np.int32(0xFF) - np.int32(round(att_dB * 8)) for att_dB in self.att_probe]
        self.att_reg = 0x00000000


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

            # sweep frequencies
            for i in range(len(self.freq_probe_scan_ftw)):
                self.core.break_realtime()

                # set freq and ampl for probe
                freq_ftw = self.freq_probe_scan_ftw[i]
                self.dds_probe.set_mu(freq_ftw, asf=self.ampl_probe_asf)
                self.core.break_realtime()

                # set att for probe
                self.dds_probe.set_att_mu(self.att_probe[i])
                self.core.break_realtime()

                # run pulse sequence from core DMA (repump on)
                self.core_dma.playback_handle(handle_on)

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, 1, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

                # run pulse sequence from core DMA (repump off)
                self.core_dma.playback_handle(handle_off)
                with parallel:
                    self.update_dataset(freq_ftw, 0, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

        # reset board profiles
        self.dds_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # sequence w/866 (cooling repump) on
        with self.core_dma.record(_DMA_HANDLE_ON):
            # turn cooling repump on
            self.dds_repump_cooling.cfg_sw(1)

            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_doppler_cooling_mu)

            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(0b0101)
                self.pmt_gating_edge(self.time_probe_mu)

        # sequence w/866 (cooling repump) off
        with self.core_dma.record(_DMA_HANDLE_OFF):
            # turn cooling repump off
            self.dds_repump_cooling.cfg_sw(0)

            # pump on, probe off
            with parallel:
                self.dds_board.cfg_switches(0b0010)
                delay_mu(self.time_doppler_cooling_mu)

            # probe on, pump off, PMT start recording
            with parallel:
                self.dds_board.cfg_switches(0b0001)
                self.pmt_gating_edge(self.time_probe_mu)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set pump waveform
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf)
        self.core.break_realtime()

        # initialize repump beam and set waveform
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf)
        self.core.break_realtime()

        # set up sampler
        self.adc.set_gain_mu(self.photodiode_channel, self.photodiode_gain)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, repump_status, pmt_counts, sampler_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('temperature_measurement', [freq_ftw * self.ftw_to_frequency, repump_status, pmt_counts, sampler_mu * self.adc_mu_to_volts])


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
            binned_count_list = np.heaviside(np.array(count_list[1]) - self.pmt_discrimination, 1)
            self.temperature_measurement_processed[i] = np.array([freq_mhz, np.mean(binned_count_list), np.std(binned_count_list)])

        print(self.temperature_measurement_processed)
