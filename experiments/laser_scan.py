import numpy as np
from artiq.experiment import *
# todo: upload data to labrad
# todo: check scannable works correctly
_DMA_HANDLE_TIMESWEEP = "timesweep_rdx"


class LaserScan(EnvExperiment):
    """
    729nm Laser Scan but better
    Gets 729nm Spectrum
    """

    #kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "dds_board_num",
        "dds_board_qubit_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_729_us",                    NumberValue(default=400, ndecimals=5, step=1, min=1, max=10000000))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(default=RangeScan(100, 120, 1001),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=3))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =              self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_cooling_mu =          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =          self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_729_mu =              self.core.seconds_to_mu(self.time_729_us * us)
        self.time_repump_qubit_mu =     self.core.seconds_to_mu(self.time_repump_qubit_us * us)

        # DDS devices
        self.dds_board =                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =          self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_pump =                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # process scan frequencies
        self.ftw_to_mhz =               1e3 / (2 ** 32 - 1)
        self.freq_qubit_scan_ftw =      [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz]

        # convert dds values to machine units - frequency
        self.freq_pump_cooling_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert dds values to machine units - amplitude
        self.ampl_pump_cooling_asf =    self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =    self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # sort out attenuation
        self.att_cooling_mu =           np.int32(0xFF) - np.int32(round(self.att_pump_cooling_dB * 8))
        self.att_readout_mu =           np.int32(0xFF) - np.int32(round(self.att_pump_readout_dB * 8))

        # set up datasets
        self.set_dataset("laser_scan_rdx", [])
        self.setattr_dataset("laser_scan_rdx")
        self.set_dataset("laser_scan_rdx_processed", np.zeros([len(self.freq_qubit_scan_ftw), 3]))
        self.setattr_dataset("laser_scan_rdx_processed")

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma
        self.DMArecord()
        self.core.break_realtime()

        # get dma handle
        handle = self.core_dma.get_handle(_DMA_HANDLE_TIMESWEEP)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # set frequencies
            for freq_ftw in self.freq_qubit_scan_ftw:
                self.core.break_realtime()

                # set waveform for qubit DDS
                self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf)
                self.core.break_realtime()

                # run sequence
                self.core_dma.playback_handle(handle)

                # update dataset
                with parallel:
                    self.update_dataset(freq_ftw, self.pmt_counter.fetch_count())
                    self.core.break_realtime()

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE_TIMESWEEP):
            # repump pulse
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

            # set cooling waveform
            self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf)

            # cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # 729
            self.dds_qubit.cfg_sw(1)
            delay_mu(self.time_729_mu)
            self.dds_qubit.cfg_sw(0)

            # set readout waveform
            self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf)

            # readout
            self.dds_board.cfg_switches(0b0110)
            self.pmt_gating_edge(self.time_readout_mu)
            self.dds_board.cfg_switches(0b0100)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set AOM DDS waveforms
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('laser_scan_rdx', [freq_ftw * self.ftw_to_mhz, pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.laser_scan_rdx = np.array(self.laser_scan_rdx)

        # get sorted x-values (frequency)
        freq_list_mhz = sorted(set(self.laser_scan_rdx[:, 0]))

        # collate results
        collated_results = {
            freq: []
            for freq in freq_list_mhz
        }
        for freq_mhz, pmt_counts in self.laser_scan_rdx:
            collated_results[freq_mhz].append(pmt_counts)

        # process counts for mean and std and put into processed dataset
        for i, (freq_mhz, count_list) in enumerate(collated_results.items()):
            self.laser_scan_rdx_processed[i] = np.array([freq_mhz, np.mean(count_list), np.std(count_list)])

        print(self.laser_scan_rdx_processed)
