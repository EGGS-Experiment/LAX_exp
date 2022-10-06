import numpy as np
from artiq.experiment import *

_DMA_HANDLE_RESET = "rabi_flopping_reset"
_DMA_HANDLE_READOUT = "rabi_flopping_readout"


class RabiFloppingRDX(EnvExperiment):
    """
    Rabi Flopping RDX
    Measures ion fluorescence vs 729nm pulse time and frequency, but better than before
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",
        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_probe_us",
        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",
        "freq_qubit_mhz",
        "ampl_probe_pct",
        "ampl_pump_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",
        "att_probe_dB",
        "att_pump_cooling_dB",
        "att_pump_readout_dB"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

        # qubit time scan
        self.setattr_argument("time_rabi_us_list",              Scannable(default=
                                                                          RangeScan(0, 200, 1001, randomize=True),
                                                                          global_min=1, global_max=100000, global_step=1,
                                                                          unit="us", scale=1, ndecimals=0
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
        self.pmt_counter =                      self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_profileswitch_delay_mu =      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_repump_qubit_mu =             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                  self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                  self.core.seconds_to_mu(self.time_readout_us * us)

        # rabi flopping timing
        self.time_rabi_mu_list =                [self.core.seconds_to_mu(time_us * us) for time_us in self.time_rabi_us_list]
        max_time_us =                           np.max(list(self.time_rabi_us_list))
        self.time_delay_mu_list =               [self.core.seconds_to_mu((max_time_us - time_us) * us) for time_us in self.time_rabi_us_list]
        self.num_time_points_list =             list(range(len(self.time_rabi_mu_list)))

        # DDS devices
        self.dds_board =                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_pump =                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_pump_cooling_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                   self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_pump_asf =                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_pct / 100)
        self.ampl_repump_cooling_asf =          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                   self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # sort out attenuation
        self.att_cooling_mu =                   np.int32(0xFF) - np.int32(round(self.att_pump_cooling_dB * 8))
        self.att_readout_mu =                   np.int32(0xFF) - np.int32(round(self.att_pump_readout_dB * 8))

        # set up datasets
        self.set_dataset("rabi_flopping_rdx", [])
        self.setattr_dataset("rabi_flopping_rdx")
        self.set_dataset("rabi_flopping_rdx_processed", np.zeros([len(self.time_rabi_mu_list), 2]))
        self.setattr_dataset("rabi_flopping_rdx_processed")

        self.set_dataset("parameters", [self.repetitions, self.freq_qubit_ftw, self.att_cooling_mu, self.att_readout_mu])

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        handle_reset = self.core_dma.get_handle(_DMA_HANDLE_RESET)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep time
            for i in self.num_time_points_list:

                # get timings
                time_delay_mu = self.time_delay_mu_list[i]
                time_rabi_mu = self.time_rabi_mu_list[i]

                # run repump and cooling
                self.core_dma.playback_handle(handle_reset)

                # wait given time
                delay_mu(time_delay_mu)

                # rabi flopping w/qubit laser
                self.dds_qubit.cfg_sw(1)
                delay_mu(time_rabi_mu)
                self.dds_qubit.cfg_sw(0)

                # do readout
                self.core_dma.playback_handle(handle_readout)

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_mu, self.pmt_counter.fetch_count())
                    self.core.break_realtime()

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE_RESET):
            with sequential:
                # qubit repump (854) pulse
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

                # set cooling waveform
                self.dds_pump.set_att_mu(self.att_cooling_mu)
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # do cooling
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
                # set readout waveform
                self.dds_pump.set_att_mu(self.att_readout_mu)
                with parallel:
                    self.dds_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # readout pulse
                self.dds_board.cfg_switches(0b0110)
                self.pmt_gating_edge(self.time_readout_mu)
                self.dds_board.cfg_switches(0b0100)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set attenuations for now so we can change them later
        att_reg_old = np.int32(self.dds_board.get_att_mu())
        self.dds_board.set_all_att_mu(att_reg_old)
        self.core.break_realtime()

        # set AOM DDS waveforms
        # profile 0 is cooling, profile 1 is readout
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()

        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.core.break_realtime()

    @rpc(flags={"async"})
    def update_dataset(self, time_mu, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('rabi_flopping_rdx', [self.core.mu_to_seconds(time_mu), pmt_counts])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.rabi_flopping_rdx = np.array(self.rabi_flopping_rdx)

        # get sorted x-values (time, seconds)
        time_list_s = sorted(set(self.rabi_flopping_rdx[:, 0]))

        # collate results
        collated_results = {
            time: []
            for time in time_list_s
        }
        for time_s, pmt_counts in self.rabi_flopping_rdx:
            collated_results[time_s].append(pmt_counts)

        # process counts for mean and std and put into processed dataset
        for i, (time_s, count_list) in enumerate(collated_results.items()):
            self.rabi_flopping_rdx_processed[i] = np.array([time_s, np.mean(count_list)])

        print(self.rabi_flopping_rdx_processed)
