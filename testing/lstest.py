import numpy as np
from artiq.experiment import *
from LAX_exp.extensions.analysis import groupBy, discriminateCounts

class analtest(EnvExperiment):
    """
    analysis test
    """

    #kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",
        "time_profileswitch_delay_us",

        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",

        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",

        "pmt_discrimination"
    ]


    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=5, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_729_us",                    NumberValue(default=40, ndecimals=5, step=1, min=1, max=10000000))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(
                                                                    default=CenterScan(104.463, 0.1, 5e-3, randomize=True),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5
                                                                ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_redist_mu =                                   self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_729_mu =                                      self.core.seconds_to_mu(self.time_729_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_profileswitch_delay_mu =                      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # process scan frequencies
        self.ftw_to_mhz =                                       1e3 / (2 ** 32 - 1)
        self.freq_qubit_scan_ftw =                              [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz]

        # convert dds values to machine units - frequency
        self.freq_redist_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert dds values to machine units - amplitude
        self.ampl_redist_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                   self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # set up datasets
        self.set_dataset("laser_scan", np.zeros([8010, 2]))
        self.setattr_dataset("laser_scan")

    def run(self):
        pass

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        laser_scan_tmp = self.tmpreplacedata()

        res_tmp = groupBy(laser_scan_tmp, combine=False)
        for i, (k, v) in enumerate(res_tmp.items()):
            self.laser_scan_processed[i, 0] = k
            self.laser_scan_processed[i, 1] = discriminateCounts(v, 17)

        from scipy.signal import find_peaks
        # todo: set height requirement, threshold, and distance
        th1 = find_peaks(self.laser_scan_processed[:, 1])[0]
        yz0 = np.zeros([len(th1),2])
        for i, val in enumerate(th1):
            yz0[i, 0] = self.laser_scan_processed[val, 0]
            yz0[i, 1] = self.laser_scan_processed[val, 1]

        print(yz0)
        self.set_dataset("laser_scan_peaks", yz0)
        self.setattr_dataset("laser_scan_peaks")

    def tmpreplacedata(self):
        import h5py
        laser_scan_tmp = np.array([])
        with h5py.File("C:\\Users\\EGGS1\\Documents\\datatmp\\000008927-LaserScan.h5", "r") as f:
            laser_scan_tmp = np.zeros([len(f['datasets']['laser_scan']), 2])
            f['datasets']['laser_scan'].read_direct(laser_scan_tmp)

        x_len = len(set(laser_scan_tmp[:, 0]))
        self.set_dataset("laser_scan_processed", np.zeros([x_len, 2]))
        self.setattr_dataset("laser_scan_processed")
        return laser_scan_tmp
