import numpy as np
from artiq.experiment import *

_DMA_HANDLE_REPUMP_ON = "temperature_measurement_repump_on"
_DMA_HANDLE_REPUMP_OFF = "temperature_measurement_repump_off"


class TemperatureMeasurement(EnvExperiment):
    """
    Temperature Measurement
    Measures ion fluorescence for a single detuning and sweeps frequency within each trial.
    Uses only the cooling beam since probe beam is polarized; keeps 866 on whole time (not once on, once off).
    """

    #kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_doppler_cooling_us",
        "time_profileswitch_delay_us",

        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "freq_pump_cooling_mhz",
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
        self.setattr_argument("repetitions",            NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_probe_us",          NumberValue(default=50, ndecimals=5, step=1, min=1, max=1000))

        # probe frequency scan
        self.setattr_argument("freq_probe_scan_mhz",    Scannable(
                                                            default=RangeScan(80, 140, 25, randomize=True),
                                                            global_min=10, global_max=200, global_step=1,
                                                            unit="MHz", scale=1, ndecimals=6
                                                        ))

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
        self.time_profileswitch_delay_mu =              self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)

        # DDS devices
        self.dds_board =                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_pump =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))

        # convert dds values to machine units - probe
        self.ftw_to_frequency =                         1e9 / (2**32 - 1)
        self.freq_probe_scan_ftw =                      [self.dds_pump.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_probe_scan_mhz]

        # convert dds values to machine units - everything else
        self.freq_pump_cooling_ftw =                    self.dds_pump.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_repump_cooling_ftw =                  self.dds_pump.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)

        self.ampl_pump_cooling_asf =                    self.dds_pump.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_repump_cooling_asf =                  self.dds_pump.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

        # ADC
        self.adc = self.get_device("sampler0")
        self.adc_buffer = np.zeros(8, dtype=int)
        self.adc_mu_to_volts = (10 ** (1 - self.photodiode_gain)) / (2 ** 15)

        # asf adjustment
        # background: 1.8uW
        asf_freqs = sorted(list(self.freq_probe_scan_mhz))
        #asf_vals = np.array([39.00, 26.60, 20.25, 16.45, 15.85, 15.00, 14.60, 14.45, 14.25, 14.00, 14.05, 14.17, 14.35, 14.93, 15.58, 16.24, 17.12, 17.95, 19.10, 20.89, 21.04, 26.95, 32.08, 39.22, 48.20])/100
        asf_vals = np.array([35,30,30,30,27.5,27.5,27.5,30,26,27.5,27.5,27.5,28,27.5,30,30,30,30,30,35,37.5,37.5,37.5,37.5,37.5])/100
        asf_dict = dict(np.concatenate([[asf_freqs], [asf_vals]]).transpose())
        self.ampl_probe_scan_asf = [self.dds_pump.amplitude_to_asf(asf_dict[freq]) for freq in self.freq_probe_scan_mhz]

        # set up datasets
        self.set_dataset("temperature_measurement", [])
        self.setattr_dataset("temperature_measurement")
        self.set_dataset("temperature_measurement_processed", np.zeros([len(self.freq_probe_scan_mhz), 3]))
        self.setattr_dataset("temperature_measurement_processed")


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
        handle_repump_on = self.core_dma.get_handle(_DMA_HANDLE_REPUMP_ON)
        handle_repump_off = self.core_dma.get_handle(_DMA_HANDLE_REPUMP_OFF)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep detuning
            for i in range(len(self.freq_probe_scan_ftw)):

                # set freq and ampl for probe
                freq_ftw = self.freq_probe_scan_ftw[i]
                ampl_asf = self.ampl_probe_scan_asf[i]
                self.dds_pump.set_mu(freq_ftw, asf=ampl_asf, profile=1)
                self.core.break_realtime()

                # todo: get photodiode reading

                # run sequence from DMA with repump on
                self.core_dma.playback_handle(handle_repump_on)

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, 1, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

                # run sequence from DMA with repump off
                self.core_dma.playback_handle(handle_repump_off)

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, 0, self.pmt_counter.fetch_count(), self.adc_buffer[self.photodiode_channel])
                    self.core.break_realtime()

        # reset board profiles & AOMs after experiment
        self.dds_board.set_profile(0)
        self.dds_board.cfg_switches(0b1110)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # sequence w/repump on
        with self.core_dma.record(_DMA_HANDLE_REPUMP_ON):
            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # doppler cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_doppler_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # read out
            self.dds_board.cfg_switches(0b0110)
            self.pmt_gating_edge(self.time_probe_mu)
            self.dds_board.cfg_switches(0b0100)


        # sequence w/repump off
        with self.core_dma.record(_DMA_HANDLE_REPUMP_OFF):
            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # doppler cooling
            self.dds_board.cfg_switches(0b0010)
            delay_mu(self.time_doppler_cooling_mu)
            self.dds_board.cfg_switches(0b0000)

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # read out
            self.dds_board.cfg_switches(0b0010)
            self.pmt_gating_edge(self.time_probe_mu)
            self.dds_board.cfg_switches(0b0000)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set pump waveform (profile 0 is doppler cooling, profile 1 is readout)
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=1)
        self.core.break_realtime()

        # initialize repump beam and set waveform
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
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
