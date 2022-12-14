import numpy as np
from random import shuffle
from artiq.experiment import *

_DMA_HANDLE_INITIALIZE = "ramsey_initialize"
_DMA_HANDLE_RAMSEY = "ramsey_pulse"
_DMA_HANDLE_READOUT = "ramsey_readout"

# todo: make max time the same as readout time
# pi/2 pulses separated by given time
# sweep frequency over the bsb


class RamseySpectroscopy(EnvExperiment):
    """
    Ramsey Spectroscopy
    Measures population transfer as a function of detuning after applying two pi/2-pulses separated by some time.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",

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
        "freq_qubit_mhz",

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
        self.setattr_argument("repetitions",                        NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # ramsey
        self.setattr_argument("time_ramsey_us",                     NumberValue(default=400, ndecimals=5, step=1, min=1, max=1000000))

        # readout
        self.setattr_argument("freq_bsb_scan_mhz",                  Scannable(
                                                                            default=
                                                                            CenterScan(105.214, 0.04, 0.001),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                    ))

        self.setattr_argument("time_half_pipulse_us",               NumberValue(default=125, ndecimals=5, step=1, min=1, max=10000))


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_doppler_cooling_mu =                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_redist_mu =                                   self.core.seconds_to_mu(self.time_redist_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_profileswitch_delay_mu =                      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_ramsey_mu =                                   self.core.seconds_to_mu(self.time_ramsey_us * us)
        self.time_half_pipulse_mu =                             self.core.seconds_to_mu(self.time_half_pipulse_us * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.ftw_to_mhz =                                       1e3 / (2 ** 32 - 1)
        self.freq_redist_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                                   self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        # process scan frequencies
        self.freq_bsb_scan_ftw =                                [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in list(self.freq_bsb_scan_mhz)]

        # convert amplitude to asf
        self.ampl_redist_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                   self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)

        # set up datasets
        self.set_dataset("ramsey_sequence", [])
        self.setattr_dataset("ramsey_sequence")
        self.set_dataset("ramsey_sequence_processed", np.zeros([len(self.freq_bsb_scan_ftw), 3]))
        self.setattr_dataset("ramsey_sequence_processed")


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
        handle_initialize = self.core_dma.get_handle(_DMA_HANDLE_INITIALIZE)
        handle_ramsey = self.core_dma.get_handle(_DMA_HANDLE_RAMSEY)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for i in range(self.repetitions):

            # sweep final pi-pulse frequency
            for freq_ftw in self.freq_bsb_scan_ftw:

                # set ramsey frequency in advance
                self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf, profile=0)
                self.core.break_realtime()

                # initialize by running doppler cooling and spin polarization
                self.core_dma.playback_handle(handle_initialize)

                # run ramsey spectroscopy
                self.core_dma.playback_handle(handle_ramsey)

                # read out
                self.core_dma.playback_handle(handle_readout)

                # record data
                self.update_dataset(freq_ftw, self.pmt_counter.fetch_count())
                self.core.break_realtime()

        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # doppler cooling sequence
        with self.core_dma.record(_DMA_HANDLE_INITIALIZE):

            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # repump qubit into S-1/2 state
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

            # doppler cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_doppler_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # spin polarization/redistribute S-1/2 (397)
            self.dds_board.cfg_switches(0b0101)
            delay_mu(self.time_redist_mu)
            self.dds_board.cfg_switches(0b0100)

        # ramsey sequence
        with self.core_dma.record(_DMA_HANDLE_RAMSEY):

            # set spectroscopy waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # first pi/2 pulse
            self.dds_qubit.cfg_sw(True)
            delay_mu(self.time_half_pipulse_mu)
            self.dds_qubit.cfg_sw(False)

            # free evolution
            delay_mu(self.time_ramsey_mu)

            # second pi/2 pulse
            self.dds_qubit.cfg_sw(True)
            delay_mu(self.time_half_pipulse_mu)
            self.dds_qubit.cfg_sw(False)

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # read out
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
        # profile 0 = cooling; profile 1 = readout (red-detuned); profile 2 = readout (blue-detuned)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=0)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=1)
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('ramsey_sequence', [freq_ftw * self.ftw_to_mhz, pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.ramsey_sequence = np.array(self.ramsey_sequence)

        # get sorted x-values (frequency)
        freq_list_mhz = sorted(set(self.ramsey_sequence[:, 0]))

        # collate results
        collated_results = {
            freq: []
            for freq in freq_list_mhz
        }
        for freq_mhz, pmt_counts in self.ramsey_sequence:
            collated_results[freq_mhz].append(pmt_counts)

        # process counts for mean and std and put into processed dataset
        for i, (freq_mhz, count_list) in enumerate(collated_results.items()):
            binned_count_list = np.heaviside(np.array(count_list) - self.pmt_discrimination, 1)
            self.ramsey_sequence_processed[i] = np.array([freq_mhz, np.mean(binned_count_list), np.std(binned_count_list)])
